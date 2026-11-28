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
              ele_anchor_pt, harmonic, fundamental_frequency, metadata=Dict{String,Any}())

A rectangular (`x, y, z`) or cylindrical (`r, theta, z`) grid of physical
electric (V/m) and magnetic (T) field samples. The grid origin, spacing,
and size follow the native axis order. Spatial coordinates are in metres and
the azimuth is in radians. The harmonic and fundamental frequency (Hz) define
the operating frequency.

Access native components as properties, for example `fm.Bz` or `fm.Etheta`.
The three-dimensional arrays are mutable. An omitted native component starts
at zero and is allocated when accessed. For construction from arrays, pass
a dictionary such as `Dict(:Bz => bz)` holding physical samples. HDF5 loading
applies the file's field scale and RF phase once.

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
    const components::Dict{Symbol, Array{<:Number, 3}}
    const metadata::Dict{String, Any}

    function FieldMesh(; geometry, grid_origin, grid_spacing, grid_size, components,
                       ele_anchor_pt, harmonic, fundamental_frequency, metadata=Dict{String,Any}())
        _check_field_components(geometry, grid_size, components)
        file_attrs = merge(Dict{String,Any}("gridLowerBound" => [0, 0, 0],
                                            "fieldScale" => 1.0, "RFphase" => 0.0), metadata)
        return new(geometry, AXIS_LABELS_FROM_GEOMETRY[geometry], grid_origin, grid_spacing,
                   grid_size, ele_anchor_pt, harmonic, fundamental_frequency, components, file_attrs)
    end
end

function _check_field_components(geometry, grid_size, components)
    native = _native_components(geometry)
    for (name, data) in components
        name in native || throw(ArgumentError("$name is not a native $geometry field component"))
        size(data) == grid_size || throw(DimensionMismatch("$name has size $(size(data)); expected $grid_size"))
    end
end

function _field_component(fm::FieldMesh, name::Symbol)
    name in _native_components(fm.geometry) ||
        throw(ArgumentError("$name is not a native $(fm.geometry) field component"))
    return get!(fm.components, name) do
        zeros(is_static(fm) ? Float64 : ComplexF64, fm.grid_size)
    end
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

function _file_field_factor(scale, rf_phase)
    return iszero(rf_phase) ? scale : scale * cispi(-2 * rf_phase)
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
    file_factor = _file_field_factor(a["fieldScale"], a["RFphase"])

    components = Dict{Symbol, Array{<:Number, 3}}()
    for (component, path) in FIELD_RECORD_FROM_COMPONENT
        haskey(h5, path) || continue
        node = h5[path]
        record = first(split(path, '/'))
        dim = Tuple(attrs(node)["unitDimension"])
        dim == EXPECTED_RECORD_UNIT_DIMENSION[record] ||
            throw(ArgumentError("$path has incompatible unitDimension $dim"))
        samples = component_data(node; axis_labels=file_axes)
        if iszero(harmonic) && !(eltype(samples) <: Real)
            samples = real.(samples)
        end
        components[component] = isone(file_factor) ? samples : file_factor .* samples
    end

    for key in ("gridLowerBound", "name", "gridCurvatureRadius", "masterParameter", "fieldScale", "RFphase")
        if haskey(a, key)
            metadata[key] = a[key]
        end
    end
    return FieldMesh(; geometry, grid_origin, grid_spacing, grid_size, ele_anchor_pt,
                     harmonic, fundamental_frequency, components, metadata)
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
    return all(!startswith(string(k), "B") || all(iszero, v) for (k, v) in fm.components)
end

"""
    is_pure_magnetic(fm::FieldMesh) -> Bool

True if there are no non-zero electric field components.
"""
function is_pure_magnetic(fm::FieldMesh)
    return all(!startswith(string(k), "E") || all(iszero, v) for (k, v) in fm.components)
end

"""
    field_magnitude(fm, field)

Pointwise magnitude of the physical electric (`:E`, V/m) or magnetic (`:B`, T)
field in the native basis. Complex RF components contribute their squared moduli.
"""
function field_magnitude(fm::FieldMesh, field::Symbol)
    field in (:E, :B) || throw(ArgumentError("Expected :E or :B, got $field"))
    # Scalar zeros broadcast over the stored arrays without allocating missing components.
    a, b, c = map(axis -> get(fm.components, Symbol(field, axis), 0), fm.axis_labels)
    if a isa Number && b isa Number && c isa Number
        return zeros(fm.grid_size)
    end
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
        "RFphase" => get(fm.metadata, "RFphase", 0.0),
        "fieldScale" => get(fm.metadata, "fieldScale", 1.0),
    ))
    # A map originally scaled to zero may now contain edited, nonzero fields.
    if iszero(attrs["fieldScale"]) && any(data -> any(!iszero, data), values(fm.components))
        attrs["fieldScale"] = 1.0
    end
    file_factor = _file_field_factor(attrs["fieldScale"], attrs["RFphase"])
    components = Dict{String, Any}()
    for (key, data) in fm.components
        path = FIELD_RECORD_FROM_COMPONENT[key]
        components[path] = isone(file_factor) || iszero(file_factor) ? data : data ./ file_factor
    end
    return Dict{String, Any}("attrs" => attrs, "components" => components)
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
    a.metadata == b.metadata || return false
    for name in _native_components(a.geometry)
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
    metadata=deepcopy(fm.metadata),
)
