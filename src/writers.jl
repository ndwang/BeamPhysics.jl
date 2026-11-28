"""
Utilities for writing particle and field data to HDF5 files in openPMD format.
"""

"""
    pmd_init(h5; basePath="/data/%T/", particlesPath="./")

Root attribute initialization for openPMD particle files.
h5 should be the root of the file.
"""
function pmd_init(h5; basePath="/data/%T/", particlesPath="./")
    attrs(h5)["basePath"] = basePath
    attrs(h5)["dataType"] = "openPMD"
    attrs(h5)["openPMD"] = "2.0.0"
    attrs(h5)["openPMDextension"] = "BeamPhysics;SpeciesType"
    attrs(h5)["particlesPath"] = particlesPath
end

"""
    pmd_field_init(h5; externalFieldPath="/ExternalFieldPath/%T/")

Root attribute initialization for an openPMD-beamphysics External Field Mesh.
h5 should be the root of the file.
"""
function pmd_field_init(h5; externalFieldPath="/ExternalFieldPath/%T/")
    attrs(h5)["dataType"] = "openPMD"
    attrs(h5)["openPMD"] = "2.0.0"
    attrs(h5)["openPMDextension"] = "BeamPhysics"
    attrs(h5)["externalFieldPath"] = externalFieldPath
end

"""
    write_component_data(h5, name, data; unit=nothing, reverse_axes=true)

Write a single record component to an HDF5 group.

If all values in `data` are equal (constant component), a group is created with
`value` and `shape` stored as attributes. Otherwise, a dataset is created.

Multidimensional datasets are written in Python axis order by default.
Set reverse_axes=false for a reversed openPMD axis layout.

Unit metadata is written via `write_unit_h5` when `unit` is provided.
"""
function write_component_data(h5, name, data; unit=nothing, reverse_axes=true)
    if !isempty(data) && all(==(first(data)), data)
        # Constant component: store as group with value/shape attributes
        g = create_group(h5, name)
        attrs(g)["value"] = first(data)
        attrs(g)["shape"] = collect(size(data))
    else
        if reverse_axes && ndims(data) > 1
            data = permutedims(data, Tuple(ndims(data):-1:1))
        end
        # HDF5 needs a dense buffer for strided or reordered selections.
        h5[name] = convert(Array, data)
        g = h5[name]
    end

    if !isnothing(unit)
        write_unit_h5(g, unit)
    end

    return g
end

function _time_offsets(t_offset, n)
    if t_offset isa Number
        return iszero(t_offset) ? nothing : fill(Float64(t_offset), n)
    end
    t_offset isa AbstractVector && length(t_offset) == n ||
        throw(DimensionMismatch("time offset must be a scalar or a vector with one value per particle"))
    return all(iszero, t_offset) ? nothing : Float64.(t_offset)
end

"""
    write_pmd_bunch(h5, data::Dict; name=nothing, t_offset=0.0)

Write particle bunch data in openPMD format.

`data` is a Dict with:
- Array keys: `"x"`, `"px"`, `"y"`, `"py"`, `"z"`, `"pz"`, `"t"`, `"status"`, `"weight"`
- String key: `"species"`
- Int key: `"n_particle"`
- Float key: `"charge"`
- Optional array key: `"id"`

Nonzero t_offset (a scalar or one value per particle, in seconds) is written as a timeOffset record.

See inverse routine: `ParticleGroup(h5)`
"""
function write_pmd_bunch(h5, data::Dict; name=nothing, t_offset=0.0)
    offsets = _time_offsets(t_offset, data["n_particle"])
    return _write_pmd_bunch(h5, data, name, offsets)
end

function _write_pmd_bunch(h5, data, name, offsets)
    g = isnothing(name) ? h5 : create_group(h5, name)

    # Create species subgroup
    species = data["species"]
    g = create_group(g, species)

    # Species attributes
    attrs(g)["speciesType"] = species
    attrs(g)["numParticles"] = data["n_particle"]
    attrs(g)["totalCharge"] = data["charge"]
    attrs(g)["chargeUnitSI"] = 1.0

    # Required components
    for key in ["x", "px", "y", "py", "z", "pz", "t", "status", "weight"]
        g2_name = COMPONENT_FROM_ALIAS[key]
        u = pg_units(key)
        write_component_data(g, g2_name, data[key]; unit=u)
    end

    # Optional id (no units)
    if haskey(data, "id")
        g["id"] = convert(Vector, data["id"])
    end

    for key in OPTIONAL_PARTICLE_COLUMNS
        alias = string(key)
        if haskey(data, alias)
            write_component_data(g, COMPONENT_FROM_ALIAS[alias], data[alias]; unit=pg_units(alias))
        end
    end

    if !isnothing(offsets)
        write_component_data(g, "timeOffset", offsets; unit=pg_units("t"))
    end
end

"""
    write_pmd_field(h5, data::Dict; name=nothing)

Write field mesh data in openPMD format.

`data` is a Dict with:
- `"attrs"`: flat dict of field attributes
- `"components"`: flat dict of component name => array

See inverse routine: `load_field_attrs`
"""
function write_pmd_field(h5, data::Dict; name=nothing)
    g = isnothing(name) ? h5 : create_group(h5, name)

    # Validate attrs
    field_attrs, other = load_field_attrs(data["attrs"])

    # Write attributes
    for (k, v) in field_attrs
        attrs(g)[k] = v
    end
    for (k, v) in other
        attrs(g)[k] = v
    end

    # Write components
    reverse_axes = !(Tuple(field_attrs["axisLabels"]) in REVERSED_FIELD_AXES)
    grid_size = Tuple(field_attrs["gridSize"])
    for (key, val) in data["components"]
        component = get(FIELD_RECORD_FROM_COMPONENT, Symbol(key), key)
        haskey(FIELD_COMPONENT_FROM_RECORD, component) ||
            throw(ArgumentError("Unknown field component: $key"))
        size(val) == grid_size ||
            throw(DimensionMismatch("$key does not match gridSize"))
        write_component_data(g, component, ComplexF64.(val);
                             unit=pg_units(String(FIELD_COMPONENT_FROM_RECORD[component])), reverse_axes)
    end
end

"""
    particle_data_dict(pg::AbstractParticleGroup) -> Dict{String, Any}

Extract a Dict from a ParticleGroup or ParticleGroupView in the format expected by `write_pmd_bunch`.
"""
function particle_data_dict(pg::AbstractParticleGroup)
    return Dict{String, Any}(
        "x" => pg.x,
        "px" => pg.px,
        "y" => pg.y,
        "py" => pg.py,
        "z" => pg.z,
        "pz" => pg.pz,
        "t" => pg.t,
        "status" => pg.status,
        "weight" => pg.weight,
        "species" => nameof(pg.species),
        "n_particle" => length(pg),
        "charge" => charge(pg),
        "id" => pg.id,
        (string(k) => pg.extra[k] for k in OPTIONAL_PARTICLE_COLUMNS if haskey(pg.extra, k))...,
    )
end

"""
    write_particle_group(filename::AbstractString, pg::AbstractParticleGroup; t_offset=0.0)

Write a ParticleGroup or ParticleGroupView to an HDF5 file in openPMD format.
This is the inverse of `ParticleGroup(filename::AbstractString)`.
"""
function write_particle_group(filename::AbstractString, pg::AbstractParticleGroup; t_offset=0.0)
    offsets = _time_offsets(t_offset, length(pg))
    data = particle_data_dict(pg)
    h5open(filename, "w") do h5
        pmd_init(h5)
        g = create_group(h5, "data/1")
        _write_pmd_bunch(g, data, nothing, offsets)
    end
end
