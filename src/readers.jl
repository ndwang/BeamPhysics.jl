"""
Utilities for reading particle and field data from HDF5 files in openPMD format.
"""

# -----------------------------------------
# Records, components, units

const PARTICLE_RECORD_COMPONENTS = Dict{String, Union{Nothing, Vector{String}}}(
    "branchIndex" => nothing,
    "chargeState" => nothing,
    "electricField" => ["x", "y", "z"],
    "elementIndex" => nothing,
    "magneticField" => ["x", "y", "z"],
    "locationInElement" => nothing,
    "momentum" => ["x", "y", "z"],
    "momentumOffset" => ["x", "y", "z"],
    "photonPolarizationAmplitude" => ["x", "y"],
    "photonPolarizationPhase" => ["x", "y"],
    "sPosition" => nothing,
    "totalMomentum" => nothing,
    "totalMomentumOffset" => nothing,
    "particleStatus" => nothing,
    "pathLength" => nothing,
    "position" => ["x", "y", "z"],
    "positionOffset" => ["x", "y", "z"],
    "spin" => ["x", "y", "z", "theta", "phi", "psi"],
    "time" => nothing,
    "timeOffset" => nothing,
    "velocity" => ["x", "y", "z"],
    "weight" => nothing
)

const FIELD_RECORD_FROM_COMPONENT = Dict(
    Symbol(prefix, axis) => record * "/" * axis
    for (prefix, record) in (("E", "electricField"), ("B", "magneticField"))
    for axis in ("x", "y", "z", "r", "theta"))
const FIELD_COMPONENT_FROM_RECORD = Dict(record => component for (component, record) in FIELD_RECORD_FROM_COMPONENT)

# Expected unit dimensions for particle and field records
const EXPECTED_RECORD_UNIT_DIMENSION = Dict{String, NTuple{7, Float64}}(
    "branchIndex" => DIMENSION["1"],
    "chargeState" => DIMENSION["1"],
    "electricField" => DIMENSION["electric_field"],
    "magneticField" => DIMENSION["magnetic_field"],
    "elementIndex" => DIMENSION["1"],
    "locationInElement" => DIMENSION["1"],
    "momentum" => DIMENSION["momentum"],
    "momentumOffset" => DIMENSION["momentum"],
    "photonPolarizationAmplitude" => DIMENSION["electric_field"],
    "photonPolarizationPhase" => DIMENSION["1"],
    "sPosition" => DIMENSION["length"],
    "totalMomentum" => DIMENSION["momentum"],
    "totalMomentumOffset" => DIMENSION["momentum"],
    "particleStatus" => DIMENSION["1"],
    "pathLength" => DIMENSION["length"],
    "position" => DIMENSION["length"],
    "positionOffset" => DIMENSION["length"],
    "spin" => DIMENSION["1"],
    "time" => DIMENSION["time"],
    "timeOffset" => DIMENSION["time"],
    "velocity" => DIMENSION["velocity"],
    "weight" => DIMENSION["charge"]
)

# Convenient aliases for components
const COMPONENT_FROM_ALIAS = Dict{String, String}(
    "t" => "time",
    "weight" => "weight",
    "status" => "particleStatus",
    "s" => "sPosition",
    "ix_branch" => "branchIndex",
    "ix_ele" => "elementIndex",
    "location" => "locationInElement",
    "p0c" => "totalMomentumOffset",
    "spin_x" => "spin/x",
    "spin_y" => "spin/y",
    "spin_z" => "spin/z",
    "charge_state" => "chargeState"
)

# Add particle position and momentum aliases.
for (g, prefix) in (("position", ""), ("momentum", "p"))
    for c in ["x", "y", "z", "r", "theta"]
        alias = prefix * c
        COMPONENT_FROM_ALIAS[alias] = g * "/" * c
    end
end

const OPTIONAL_PARTICLE_COLUMNS = (:s, :ix_branch, :ix_ele, :location, :p0c,
                                   :spin_x, :spin_y, :spin_z, :charge_state)
const REVERSED_FIELD_AXES = (("z", "y", "x"), ("z", "theta", "r"))

# Some HDF5 writers encode scalar attributes as one-element arrays.
_scalar_attr(x::AbstractArray) = only(x)
_scalar_attr(x) = x

"""
    particle_paths(h5; key="particlesPath")

Uses the basePath and particlesPath to find where openPMD particles should be.
"""
function particle_paths(h5::Union{HDF5.File,HDF5.Group}; key="particlesPath")
    basePath = String(_scalar_attr(attrs(h5)["basePath"]))
    particlesPath = String(_scalar_attr(attrs(h5)[key]))

    # "." and "./" mean current directory in openPMD; HDF5.jl cannot resolve them
    if particlesPath in (".", "./")
        particlesPath = ""
    end

    if !occursin("%T", basePath)
        return [basePath * particlesPath]
    end

    path1, path2 = split(basePath, "%T")
    prefix = strip(path1, '/')
    tlist = isempty(prefix) ? keys(h5) : keys(h5[prefix])
    return [path1 * t * path2 * particlesPath for t in tlist]
end

"""
    field_paths(h5; key="externalFieldPath")

Looks for the External Fields.
"""
function field_paths(h5; key="externalFieldPath")
    if !haskey(attrs(h5), key)
        return String[]
    end

    fpath = String(_scalar_attr(attrs(h5)[key]))

    if !occursin("%T", fpath)
        return [fpath]
    end

    path1, path2 = split(fpath, "%T")
    tlist = keys(h5[path1])
    return [path1 * t * path2 for t in tlist]
end

"""
    is_constant_component(h5)

Constant record component should have 'value' and 'shape'.
"""
function is_constant_component(h5)
    return haskey(attrs(h5), "value") && haskey(attrs(h5), "shape")
end

"""
    constant_component_value(h5)

Get value from constant component.
"""
function constant_component_value(h5)
    unitSI = _scalar_attr(attrs(h5)["unitSI"])
    val = _scalar_attr(attrs(h5)["value"])
    return unitSI == 1.0 ? val : val * unitSI
end

"""
    component_data(h5; slice=:, unit_factor=1, axis_labels=nothing)

Returns an array from an h5 component.

Parameters:
- h5: HDF5.Dataset or HDF5.Group
- slice: Slice or tuple of slices to retrieve parts of the array
- unit_factor: Additional factor to convert from SI units to output units
- axis_labels: Required for multidimensional arrays. Supported options:
  * ("z", "y", "x")
  * ("z", "theta", "r")
  * ("x", "y", "z")
  * ("r", "theta", "z")
"""
function component_data(h5; slice=:, unit_factor=1, axis_labels=nothing)
    factor = Float64(_scalar_attr(get(attrs(h5), "unitSI", 1))) * unit_factor
    if is_constant_component(h5)
        dims = Tuple(Int.(attrs(h5)["shape"]))
        indices = slice isa Tuple ? slice : (slice,)
        indices = (indices..., ntuple(_ -> Colon(), length(dims) - length(indices))...)
        checkbounds(CartesianIndices(dims), indices...)
        shape = Tuple(i isa Colon ? dims[k] : length(i)
                      for (k, i) in enumerate(indices) if !(i isa Integer))
        value = _scalar_attr(attrs(h5)["value"])
        dat = isempty(shape) ? value : fill(value, shape)
    elseif ndims(h5) == 1
        dat = slice isa Colon ? read(h5) : h5[slice]
    else
        isnothing(axis_labels) && throw(ArgumentError("axis_labels required for multidimensional arrays"))
        # HDF5.jl exposes the file axes in reverse. A reversed openPMD field
        # layout cancels that reversal; all other layouts need one transpose.
        legacy_f = _scalar_attr(get(attrs(h5), "gridDataOrder", "")) == "F"
        reverse_axes = !(legacy_f || Tuple(axis_labels) in REVERSED_FIELD_AXES)
        if slice isa Colon
            dat = read(h5)
        else
            indices = slice isa Tuple ? slice : (slice,)
            indices = (indices..., ntuple(_ -> Colon(), ndims(h5) - length(indices))...)
            dat = h5[(reverse_axes ? reverse(indices) : indices)...]
        end
        if reverse_axes && dat isa AbstractArray && ndims(dat) > 1
            dat = permutedims(dat, Tuple(ndims(dat):-1:1))
        end
    end

    factor == 1 && return dat
    T = eltype(dat)
    return dat .* (T <: Union{AbstractFloat,Complex} ? convert(T, factor) : factor)
end

"""
    offset_component_name(component_name)

Many components can also have an offset, as in:
    position/x
    positionOffset/c
Return the appropriate name.
"""
function offset_component_name(component_name)
    x = split(component_name, "/")
    if length(x) == 1
        return x[1] * "Offset"
    else
        return x[1] * "Offset/" * x[2]
    end
end

"""
    particle_array(h5, component; slice=:, include_offset=true)

Main routine to return particle arrays in fixed units.
Momentum is returned in eV/c and spin in units of hbar. Other units are SI.
Floating-point precision is preserved when converting units.
"""
function particle_array(h5, component; slice=:, include_offset=true)
    # Handle aliases
    if haskey(COMPONENT_FROM_ALIAS, component)
        component = COMPONENT_FROM_ALIAS[component]
    end

    if component in ("momentum/x", "momentum/y", "momentum/z", "totalMomentum", "totalMomentumOffset")
        unit_factor = C_LIGHT * EV_PER_J  # convert J/(m/s) to eV/c
    elseif startswith(component, "spin/") &&
           Tuple(get(attrs(h5[component]), "unitDimension", DIMENSION["1"])) == (2, 1, -1, 0, 0, 0, 0)
        unit_factor = inv(H_BAR * J_PER_EV)
    else
        unit_factor = 1
    end

    # Get data
    dat = component_data(h5[component], slice=slice, unit_factor=unit_factor)

    # Look for offset component
    ocomponent = offset_component_name(component)
    if include_offset && haskey(h5, ocomponent)
        offset = component_data(h5[ocomponent], slice=slice, unit_factor=unit_factor)
        dat = dat .+ offset
    end

    return dat
end

"""
    all_components(h5)

Look for possible components in a particle group.
"""
function all_components(h5)
    components = Vector{String}()
    for record_name in keys(h5)
        if !haskey(PARTICLE_RECORD_COMPONENTS, record_name)
            continue
        end

        # Look for components
        possible_components = PARTICLE_RECORD_COMPONENTS[record_name]

        if isnothing(possible_components)
            # Record is a component
            push!(components, record_name)
        else
            g = h5[record_name]
            for cname in possible_components
                if haskey(g, cname)
                    push!(components, record_name * "/" * cname)
                end
            end
        end
    end
    return components
end

"""
    component_str(particle_group, name)

Informational string from a component in a particle group (h5).
"""
function component_str(particle_group, name)
    g = particle_group[name]
    record_name = split(name, "/")[1]
    expected_dimension = EXPECTED_RECORD_UNIT_DIMENSION[record_name]
    this_dimension = NTuple{7,Float64}(attrs(g)["unitDimension"])
    dname = dimension_name(this_dimension)
    symbol = SI_SYMBOL[dname]

    s = name * " "

    if is_constant_component(g)
        val = constant_component_value(g)
        shape = attrs(g)["shape"]
        s *= "[constant $val with shape $shape]"
    else
        s *= "[" * string(length(g)) * " items]"
    end

    if symbol != "1"
        s *= " is a $dname with units: $symbol"
    end

    if expected_dimension != this_dimension
        s *= ", but expected units: " * SI_SYMBOL[dimension_name(expected_dimension)]
    end

    return s
end

# ----------------------------------
# Fields

const REQUIRED_FIELD_ATTRS = [
    # strings
    "eleAnchorPt",
    "gridGeometry",
    "axisLabels",
    # reals and ints
    "gridLowerBound",
    "gridOriginOffset",
    "gridSpacing",
    "gridSize",
    "harmonic"
]

# Dict with options
const OPTIONAL_FIELD_ATTRS = Dict{String, Any}(
    "name" => nothing,
    "gridCurvatureRadius" => nothing,
    "fundamentalFrequency" => 0,
    "RFphase" => 0,
    "fieldScale" => 1.0,
    "masterParameter" => nothing
)

"""
    load_field_attrs(attr; verbose=false)

Loads FieldMesh required and optional attributes from a dict_like object.
Non-standard attributes will be collected in an 'other' dict.

Returns tuple: (attrs, other)
"""
function load_field_attrs(attr; verbose=false)
    # Get all attrs. Will pop.
    a = Dict(attr)

    attrs = Dict{String, Any}()
    other = Dict{String, Any}()

    # Required
    for k in REQUIRED_FIELD_ATTRS
        value = pop!(a, k)
        attrs[k] = k in ("eleAnchorPt", "gridGeometry", "harmonic") ? _scalar_attr(value) : value
    end

    # Optional, filling in some defaults
    for (k, v) in OPTIONAL_FIELD_ATTRS
        if haskey(a, k)
            attrs[k] = _scalar_attr(pop!(a, k))
        elseif !isnothing(v)
            attrs[k] = v
        end
    end

    # Collect other
    for (k, v) in a
        other[k] = v
        if verbose
            println("Nonstandard attr: ", k, " ", v)
        end
    end

    return attrs, other
end
