# ── Geometry helpers ──────────────────────────────────────────────────────

const AXIS_LABELS_FROM_GEOMETRY = Dict(
    "rectangular" => (:x, :y, :z),
    "cylindrical" => (:r, :theta, :z),
)

function _native_components(geometry)
    geometry == "rectangular" && return (:Ex, :Ey, :Ez, :Bx, :By, :Bz)
    geometry == "cylindrical" && return (:Er, :Etheta, :Ez, :Br, :Btheta, :Bz)
    throw(ArgumentError("Unknown field geometry: $geometry"))
end

"""
    FieldMesh(filename::AbstractString)
    FieldMesh(group::HDF5.Group)
    FieldMesh(; geometry, grid_origin, grid_spacing, grid_size, components,
              ele_anchor_pt, harmonic, fundamental_frequency, field_scale=1.0,
              rf_phase=0.0, component_units=Dict(), metadata=Dict())

A rectangular (`x, y, z`) or cylindrical (`r, theta, z`) grid of physical
electric (V/m) and magnetic (T) field samples. The grid origin, spacing,
and size follow the native axis order. Spatial coordinates are in metres and
the azimuth is in radians. The harmonic and fundamental frequency (Hz) define
the operating frequency.

Access native components as properties, for example `fm.Bz` or `fm.Etheta`.
These are writable physical-field views of the unscaled arrays in `fm.components`.
Set `fm.field_scale` to change the overall amplitude and `fm.rf_phase` to
change phase in turns. Edit `fm.components[:Bz]` to change the underlying profile,
including while the scale is zero. Omitted native components read as zero.
Constructor arrays use SI units by default. Use `component_units` to specify
the units of stored samples.

```julia
fm = FieldMesh("solenoid.h5")
fm.Bz               # physical longitudinal magnetic field [T]
coord_vec(fm, :z)   # grid coordinates [m]
```
"""
mutable struct FieldMesh
    const geometry::String
    const axis_labels::NTuple{3, Symbol}
    grid_origin::NTuple{3, Float64}
    grid_spacing::NTuple{3, Float64}
    const grid_size::NTuple{3, Int}
    ele_anchor_pt::String
    harmonic::Int
    fundamental_frequency::Float64
    field_scale::Float64
    rf_phase::Float64
    const components::Dict{Symbol, Array{<:Number, 3}}
    const component_units::Dict{Symbol, PMDUnit}
    const metadata::Dict{String, Any}

    function FieldMesh(; geometry, grid_origin, grid_spacing, grid_size, components,
                       ele_anchor_pt, harmonic, fundamental_frequency, metadata=Dict{String,Any}(),
                       field_scale=get(metadata, "fieldScale", 1.0), rf_phase=get(metadata, "RFphase", 0.0),
                       component_units=Dict{Symbol,PMDUnit}())
        _check_field_components(geometry, grid_size, components)
        profiles = Dict{Symbol, Array{<:Number,3}}()
        for (name, data) in components
            # Physical writes can be fractional, and RF writes can change either quadrature.
            T = float(eltype(data))
            sample_type = iszero(harmonic) ? T : complex(T)
            profiles[name] = convert(Array{sample_type,3}, data)
        end
        for (name, unit) in component_units
            name in _native_components(geometry) || throw(ArgumentError("Unknown $geometry component $name"))
            unit.dims == pg_units(String(name)).dims ||
                throw(ArgumentError("$name has incompatible unitDimension $(unit.dims)"))
        end
        file_attrs = merge(Dict{String,Any}("gridLowerBound" => [0, 0, 0]), metadata)
        delete!(file_attrs, "fieldScale")
        delete!(file_attrs, "RFphase")
        return new(geometry, AXIS_LABELS_FROM_GEOMETRY[geometry], grid_origin, grid_spacing,
                   grid_size, ele_anchor_pt, harmonic, fundamental_frequency, field_scale, rf_phase,
                   profiles, component_units, file_attrs)
    end
end

function _check_field_components(geometry, grid_size, components)
    native = _native_components(geometry)
    for (name, data) in components
        name in native || throw(ArgumentError("$name is not a native $geometry field component"))
        size(data) == grid_size || throw(DimensionMismatch("$name has size $(size(data)); expected $grid_size"))
    end
end

# Look up scale and phase on the mesh so a retained view follows a voltage or phase scan.
struct PhysicalFieldView{T,A} <: AbstractArray{T,3}
    mesh::FieldMesh
    component::Symbol
    samples::A
    unit_scale::Float64
end

Base.size(v::PhysicalFieldView) = v.mesh.grid_size
Base.IndexStyle(::Type{<:PhysicalFieldView}) = IndexLinear()
_samples(v::PhysicalFieldView) = v.samples
_samples(v::PhysicalFieldView{T,Nothing}) where {T} = get(v.mesh.components, v.component, nothing)
function _field_gain(v::PhysicalFieldView{T}) where {T}
    scale = v.mesh.field_scale * v.unit_scale
    phase = v.mesh.rf_phase
    iszero(phase) && return scale
    # A retained real view of a static map continues to give the in-phase field.
    return T <: Real ? scale * cospi(2 * phase) : scale * cispi(-2 * phase)
end

@inline function Base.getindex(v::PhysicalFieldView{T}, i::Int) where {T}
    @boundscheck checkbounds(v, i)
    samples = _samples(v)
    isnothing(samples) && return zero(eltype(v))
    value = @inbounds samples[i]
    profile = is_static(v.mesh) ? real(value) : value
    return convert(T, _field_gain(v) * profile)
end

function Base.setindex!(v::PhysicalFieldView, value, i::Int)
    @boundscheck checkbounds(v, i)
    gain = _field_gain(v)
    if iszero(gain)
        iszero(value) && return value
        throw(ArgumentError("Set a nonzero field scale before writing a nonzero physical field, or edit fm.components directly"))
    end
    samples = _samples(v)
    if isnothing(samples)
        iszero(value) && return value
        samples = v.mesh.components[v.component] =
            zeros(is_static(v.mesh) ? Float64 : ComplexF64, size(v))
    end
    sample = value / gain
    if is_static(v.mesh) && eltype(samples) <: Complex
        # Static files can store a complex record. Only its real channel is physical.
        sample = complex(convert(real(eltype(samples)), sample), imag(samples[i]))
    end
    samples[i] = sample
    return value
end

function _field_component(fm::FieldMesh, name::Symbol)
    name in _native_components(fm.geometry) ||
        throw(ArgumentError("$name is not a native $(fm.geometry) field component"))
    samples = get(fm.components, name, nothing)
    unit = get(fm.component_units, name, nothing)
    T = is_static(fm) && iszero(fm.rf_phase) ? Float64 : ComplexF64
    return PhysicalFieldView{T,typeof(samples)}(fm, name, samples, isnothing(unit) ? 1.0 : unit.scale)
end

function Base.getproperty(fm::FieldMesh, name::Symbol)
    hasfield(FieldMesh, name) && return getfield(fm, name)
    return _field_component(fm, name)
end

Base.propertynames(fm::FieldMesh, private::Bool=false) =
    (fieldnames(FieldMesh)..., _native_components(fm.geometry)...)

# ── Constructors ──────────────────────────────────────────────────────────

"""
    FieldMesh(filename::AbstractString)

Load a FieldMesh from an openPMD-beamphysics HDF5 file.
"""
function FieldMesh(filename::AbstractString)
    h5open(filename, "r") do h5
        fp = field_paths(h5)
        isempty(fp) && error("No field paths found in $filename")
        length(fp) > 1 && error("Multiple field paths in $filename: $fp. Use FieldMesh(h5group) to select one.")
        return FieldMesh(h5[fp[1]])
    end
end

"""
    FieldMesh(h5::HDF5.Group)

Load a FieldMesh from an open HDF5 group.
"""
function FieldMesh(h5::HDF5.Group)
    a, metadata = load_field_attrs(attrs(h5))

    geometry = String(a["gridGeometry"])
    file_axes = NTuple{3, String}(String.(a["axisLabels"]))
    grid_origin = NTuple{3, Float64}(a["gridOriginOffset"])
    grid_spacing = NTuple{3, Float64}(a["gridSpacing"])
    grid_size = NTuple{3, Int}(a["gridSize"])
    ele_anchor_pt = String(a["eleAnchorPt"])
    harmonic = Int(a["harmonic"])
    fundamental_frequency = Float64(a["fundamentalFrequency"])
    component_units = Dict{Symbol, PMDUnit}()
    components = Dict{Symbol, Array{<:Number, 3}}()
    for (component, path) in FIELD_RECORD_FROM_COMPONENT
        haskey(h5, path) || continue
        node = h5[path]
        component_units[component] = read_unit_h5(node)
        components[component] = _component_samples(node, :, file_axes)
    end

    for key in ("gridLowerBound", "name", "gridCurvatureRadius", "masterParameter")
        if haskey(a, key)
            metadata[key] = a[key]
        end
    end
    return FieldMesh(; geometry, grid_origin, grid_spacing, grid_size, ele_anchor_pt,
                     harmonic, fundamental_frequency, field_scale=a["fieldScale"],
                     rf_phase=a["RFphase"], components, component_units, metadata)
end

# ── Grid accessors ────────────────────────────────────────────────────────

"""
    axis_index(fm::FieldMesh, key::Symbol) -> Int

Return the axis index (1-based) for a named axis label.
"""
function axis_index(fm::FieldMesh, key::Symbol)
    for i in 1:3
        fm.axis_labels[i] == key && return i
    end
    throw(ArgumentError("Axis not found: $key. Available: $(fm.axis_labels)"))
end

function _axis_range(fm::FieldMesh, i::Int)
    return range(fm.grid_origin[i]; step=fm.grid_spacing[i], length=fm.grid_size[i])
end

"""
    mins(fm::FieldMesh) -> NTuple{3, Float64}

Grid minimum coordinates.
"""
mins(fm::FieldMesh) = ntuple(i -> first(_axis_range(fm, i)), 3)

"""
    deltas(fm::FieldMesh) -> NTuple{3, Float64}

Grid spacing.
"""
deltas(fm::FieldMesh) = fm.grid_spacing

"""
    maxs(fm::FieldMesh) -> NTuple{3, Float64}

Grid maximum coordinates.
"""
maxs(fm::FieldMesh) = ntuple(i -> last(_axis_range(fm, i)), 3)

"""
    coord_vec(fm::FieldMesh, label::Symbol) -> Vector{Float64}

Coordinate vector along the named axis.
"""
function coord_vec(fm::FieldMesh, label::Symbol)
    return collect(_axis_range(fm, axis_index(fm, label)))
end

"""
    coord_vecs(fm::FieldMesh) -> Tuple of Vector{Float64}

Coordinate vectors for all three axes.
"""
coord_vecs(fm::FieldMesh) = Tuple(coord_vec(fm, l) for l in fm.axis_labels)

# Named shortcuts select from the same grid quantities as the general accessors.
for axis in (:x, :y, :z, :r, :theta)
    for (name, accessor, description) in (
        (Symbol(:d, axis), :deltas, "spacing"),
        (Symbol(axis, :min), :mins, "minimum"),
        (Symbol(axis, :max), :maxs, "maximum"),
    )
        doc = "    $name(fm::FieldMesh)\n\nGrid $description along the $axis axis."
        @eval @doc $doc $name(fm::FieldMesh) = $accessor(fm)[axis_index(fm, $(QuoteNode(axis)))]
    end
end

# ── Min/Max setters ───────────────────────────────────────────────────────

"""
    set_min!(fm::FieldMesh, label::Symbol, val::Real)

Set the grid minimum along `label` by shifting the grid origin.
"""
function set_min!(fm::FieldMesh, label::Symbol, val::Real)
    i = axis_index(fm, label)
    origin = collect(fm.grid_origin)
    origin[i] = Float64(val)
    fm.grid_origin = NTuple{3, Float64}(origin)
    return fm
end

"""
    set_max!(fm::FieldMesh, label::Symbol, val::Real)

Set the grid maximum along `label` by shifting the grid origin.
"""
function set_max!(fm::FieldMesh, label::Symbol, val::Real)
    i = axis_index(fm, label)
    current_max = maxs(fm)[i]
    origin = collect(fm.grid_origin)
    origin[i] += Float64(val) - current_max
    fm.grid_origin = NTuple{3, Float64}(origin)
    return fm
end

"""
    frequency(fm::FieldMesh) -> Float64

Oscillation frequency in Hz. Returns 0 for static fields.
"""
frequency(fm::FieldMesh) = is_static(fm) ? 0.0 : fm.harmonic * fm.fundamental_frequency

# ── Boolean queries ───────────────────────────────────────────────────────

"""
    is_static(fm::FieldMesh) -> Bool

True if the field is static (harmonic == 0).
"""
is_static(fm::FieldMesh) = fm.harmonic == 0

"""
    is_pure_electric(fm::FieldMesh) -> Bool

True if there are no non-zero magnetic field components.
"""
function is_pure_electric(fm::FieldMesh)
    return all(k -> !startswith(string(k), "B") || all(iszero, _field_component(fm, k)), keys(fm.components))
end

"""
    is_pure_magnetic(fm::FieldMesh) -> Bool

True if there are no non-zero electric field components.
"""
function is_pure_magnetic(fm::FieldMesh)
    return all(k -> !startswith(string(k), "E") || all(iszero, _field_component(fm, k)), keys(fm.components))
end

"""
    field_magnitude(fm, field)

Pointwise magnitude of the physical electric (`:E`, V/m) or magnetic (`:B`, T)
field in the native basis. Complex RF components contribute their squared moduli.
"""
function field_magnitude(fm::FieldMesh, field::Symbol)
    field in (:E, :B) || throw(ArgumentError("Expected :E or :B, got $field"))
    a, b, c = map(axis -> _field_component(fm, Symbol(field, axis)), fm.axis_labels)
    return hypot.(a, b, c)
end

"""
    interpolator(fm::FieldMesh, component::Symbol)
    interpolator(fm::FieldMesh, data::AbstractArray{<:Number,3})

Build a reusable linear interpolator for a physical field component or an array
on the same grid, such as `abs.(fm.Ez)`. Call it as `itp(x, y, z)`
for rectangular meshes or `itp(r, theta, z)` for cylindrical meshes. A
singleton axis accepts only its stored coordinate. Queries outside the grid
throw `BoundsError`.

Rebuild the interpolator after changing the mesh.
"""
interpolator(fm::FieldMesh, component::Symbol) = interpolator(fm, _field_component(fm, component))

function interpolator(fm::FieldMesh, data::AbstractArray{<:Number,3})
    size(data) == fm.grid_size || throw(DimensionMismatch("field array must match grid size $(fm.grid_size)"))
    ranges = ntuple(i -> _axis_range(fm, i), 3)
    singleton = map(r -> length(r) == 1, ranges)
    modes = map(s -> s ? Interpolations.NoInterp() :
                         Interpolations.BSpline(Interpolations.Linear()), singleton)
    grid = ntuple(3) do i
        singleton[i] ? axes(data, i) : ranges[i]
    end
    itp = Interpolations.scale(Interpolations.interpolate(data, modes), grid...)
    return function (x::Real, y::Real, z::Real)
        point = (x, y, z)
        indices = ntuple(3) do i
            singleton[i] || return point[i]
            point[i] == first(ranges[i]) || throw(BoundsError(itp, point))
            return 1
        end
        return itp(indices...)
    end
end

"""
    interpolate(fm::FieldMesh, component_or_array, point::NTuple{3,Real})

Linearly interpolate a field at one point in the mesh's native coordinates.
Use `interpolator` to evaluate multiple points without rebuilding.
"""
interpolate(fm::FieldMesh, component_or_array, point::NTuple{3,Real}) =
    interpolator(fm, component_or_array)(point...)

# ── I/O ───────────────────────────────────────────────────────────────────

"""
    _to_data_dict(fm::FieldMesh) -> Dict

Convert FieldMesh to the data dict format expected by `write_pmd_field`.
"""
function _to_data_dict(fm::FieldMesh)
    _check_field_components(fm.geometry, fm.grid_size, fm.components)
    attrs = merge(fm.metadata, Dict{String, Any}(
        "gridGeometry" => fm.geometry,
        "axisLabels" => String.(collect(fm.axis_labels)),
        "gridOriginOffset" => collect(fm.grid_origin),
        "gridSpacing" => collect(fm.grid_spacing),
        "gridSize" => collect(fm.grid_size),
        "gridLowerBound" => get(fm.metadata, "gridLowerBound", [0, 0, 0]),
        "eleAnchorPt" => fm.ele_anchor_pt,
        "harmonic" => fm.harmonic,
        "fundamentalFrequency" => fm.fundamental_frequency,
        "RFphase" => fm.rf_phase,
        "fieldScale" => fm.field_scale,
    ))
    return Dict{String, Any}("attrs" => attrs, "components" => fm.components,
                             "component_units" => fm.component_units)
end

"""
    write_fieldmesh(filename::AbstractString, fm::FieldMesh)

Write a FieldMesh to a new openPMD-beamphysics HDF5 file.
"""
function write_fieldmesh(filename::AbstractString, fm::FieldMesh)
    data = _to_data_dict(fm)
    h5open(filename, "w") do h5
        pmd_field_init(h5)
        g = create_group(h5, "ExternalFieldPath/1")
        write_pmd_field(g, data)
    end
end

"""
    write_fieldmesh(h5, fm::FieldMesh; name=nothing)

Write a FieldMesh to an open HDF5 group.
"""
function write_fieldmesh(h5, fm::FieldMesh; name=nothing)
    write_pmd_field(h5, _to_data_dict(fm); name=name)
end

# ── Base overloads ────────────────────────────────────────────────────────

function Base.show(io::IO, fm::FieldMesh)
    print(io, "FieldMesh($(fm.geometry), $(fm.grid_size))")
end

function Base.:(==)(a::FieldMesh, b::FieldMesh)
    a.geometry == b.geometry || return false
    a.axis_labels == b.axis_labels || return false
    a.grid_origin == b.grid_origin || return false
    a.grid_spacing == b.grid_spacing || return false
    a.grid_size == b.grid_size || return false
    a.ele_anchor_pt == b.ele_anchor_pt || return false
    a.harmonic == b.harmonic || return false
    a.fundamental_frequency == b.fundamental_frequency || return false
    a.field_scale == b.field_scale && a.rf_phase == b.rf_phase || return false
    a.metadata == b.metadata || return false
    for name in _native_components(a.geometry)
        default_unit = pg_units(String(name))
        get(a.component_units, name, default_unit) == get(b.component_units, name, default_unit) || return false
        x = get(a.components, name, nothing)
        y = get(b.components, name, nothing)
        x === y && continue
        if isnothing(x) || isnothing(y)
            all(iszero, something(x, y)) || return false
        elseif x != y
            return false
        end
    end
    return true
end

Base.copy(fm::FieldMesh) = FieldMesh(;
    fm.geometry, fm.grid_origin, fm.grid_spacing, fm.grid_size, fm.ele_anchor_pt,
    fm.harmonic, fm.fundamental_frequency, components=deepcopy(fm.components),
    fm.field_scale, fm.rf_phase, component_units=copy(fm.component_units), metadata=deepcopy(fm.metadata),
)
