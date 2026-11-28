abstract type AbstractParticleGroup{T<:Real, I<:Integer} end

"""
    ParticleGroup{T, I, AT, AI, AId}

Stores the particle columns of a single-species beam.
Particle columns are stored as `AbstractVector` subtypes.

# Fields
- `x,y,z::AT`: positions in [m]
- `px,py,pz::AT`: momenta in [eV/c]
- `t::AT`: time in [s]
- `status::AI`: particle status (1 for alive, 0 for dead)
- `weight::AT`: macro-particle charge in [C]
- `species::Species`: particle species
- `id::AId`: particle IDs
- `extra`: optional per-particle openPMD records
"""

struct ParticleGroup{T<:Real, I<:Integer,
                     AT<:AbstractVector{T},
                     AI<:AbstractVector{I},
                     AId<:AbstractVector{<:Integer}} <: AbstractParticleGroup{T,I}
    x::AT
    px::AT
    y::AT
    py::AT
    z::AT
    pz::AT
    t::AT
    status::AI
    weight::AT
    species::Species
    id::AId
    extra::Dict{Symbol, AbstractVector}

    function ParticleGroup{T,I,AT,AI,AId}(x, px, y, py, z, pz, t, status, weight, species, id, extra) where {T,I,AT,AI,AId}
        pg = new{T,I,AT,AI,AId}(x, px, y, py, z, pz, t, status, weight, species, id, extra)
        n = length(pg.x)
        for name in (:px, :y, :py, :z, :pz, :t, :status, :weight, :id)
            len = length(getfield(pg, name))
            len == n || throw(DimensionMismatch("column $name has length $len; expected $n to match x"))
        end
        for (name, column) in pg.extra
            len = length(column)
            len == n || throw(DimensionMismatch("extra column $name has length $len; expected $n to match x"))
        end
        return pg
    end
end

function ParticleGroup(x::AT, px::AT, y::AT, py::AT, z::AT, pz::AT, t::AT,
                       status::AI, weight::AT, species::Species, id::AId,
                       extra::Dict{Symbol, AbstractVector}) where {
                       T<:Real, I<:Integer, AT<:AbstractVector{T},
                       AI<:AbstractVector{I}, AId<:AbstractVector{<:Integer}}
    return ParticleGroup{T,I,AT,AI,AId}(x, px, y, py, z, pz, t, status, weight, species, id, extra)
end

"""
    ParticleGroup(; x, px, y, py, z, pz, t, status, weight, species, id=collect(1:length(x)), extra=Dict{Symbol, AbstractVector}())

Construct a particle group from named columns. Every column must have one entry
per particle, including IDs and optional columns. Arrays are stored without copying.
"""
function ParticleGroup(; x, px, y, py, z, pz, t, status, weight, species,
                       id=collect(1:length(x)), extra=Dict{Symbol, AbstractVector}())
    return ParticleGroup(x, px, y, py, z, pz, t, status, weight, species, id,
                         Dict{Symbol, AbstractVector}(extra))
end
"""
    ParticleGroupView

Non-owning view into a subset of particles from a `ParticleGroup` (or another
`AbstractParticleGroup`). All array fields are `SubArray`s referencing the
parent storage, so mutating a `ParticleGroupView` mutates the original data.
"""
struct ParticleGroupView{T<:Real, I<:Integer,
                         AT<:AbstractVector{T},
                         AI<:AbstractVector{I},
                         AId<:AbstractVector{<:Integer}} <: AbstractParticleGroup{T,I}
    x::AT
    px::AT
    y::AT
    py::AT
    z::AT
    pz::AT
    t::AT
    status::AI
    weight::AT
    species::Species
    id::AId
    extra::Dict{Symbol, AbstractVector}
end

function ParticleGroupView(pg::AbstractParticleGroup{T,I}, Iinds) where {T<:Real, I<:Integer}
    return ParticleGroupView(
        @view(pg.x[Iinds]), @view(pg.px[Iinds]),
        @view(pg.y[Iinds]), @view(pg.py[Iinds]),
        @view(pg.z[Iinds]), @view(pg.pz[Iinds]),
        @view(pg.t[Iinds]), @view(pg.status[Iinds]),
        @view(pg.weight[Iinds]), pg.species,
        @view(pg.id[Iinds]), Dict{Symbol, AbstractVector}(k => view(v, Iinds) for (k, v) in pg.extra),
    )
end

Base.view(pg::AbstractParticleGroup, Iinds) = ParticleGroupView(pg, Iinds)
Base.view(pg::AbstractParticleGroup, i::Integer) = view(pg, i:i)

"""
    single_particle(; x=0, px=0, y=0, py=0, z=0, pz=0, t=0, weight=1, status=1, species=Species("electron"))

Construct a `ParticleGroup` containing a single particle. All keyword arguments
are optional and default to zero (or 1 for `weight` and `status`).

# Example
```julia
pg = single_particle(z=0.01, pz=1e6)
```
"""
function single_particle(;
    x=0.0, px=0.0, y=0.0, py=0.0, z=0.0, pz=0.0,
    t=0.0, weight=1.0, status=1, species=Species("electron")
)
    x, px, y, py, z, pz, t, weight = promote(x, px, y, py, z, pz, t, weight)
    return ParticleGroup(
        [x], [px], [y], [py], [z], [pz], [t],
        [status], [weight], species, [1], Dict{Symbol, AbstractVector}()
    )
end

# Base methods
Base.length(pg::AbstractParticleGroup) = length(pg.x)

# ---- Standalone derived-property functions ----

nalive(pg::AbstractParticleGroup) = count(pg.status .== 1)
ndead(pg::AbstractParticleGroup) = length(pg) - nalive(pg)
mass(pg::AbstractParticleGroup) = massof(pg.species)
species_charge(pg::AbstractParticleGroup) = chargeof(pg.species; C=true)
charge(pg::AbstractParticleGroup) = sum(pg.weight)

momentum(pg::AbstractParticleGroup) = sqrt.(pg.px.^2 .+ pg.py.^2 .+ pg.pz.^2)
energy(pg::AbstractParticleGroup) = sqrt.(pg.px.^2 .+ pg.py.^2 .+ pg.pz.^2 .+ mass(pg)^2)
kinetic_energy(pg::AbstractParticleGroup) = energy(pg) .- mass(pg)

xp(pg::AbstractParticleGroup) = pg.px ./ pg.pz
yp(pg::AbstractParticleGroup) = pg.py ./ pg.pz

r(pg::AbstractParticleGroup) = hypot.(pg.x, pg.y)
theta(pg::AbstractParticleGroup) = atan.(pg.y, pg.x)

function pr(pg::AbstractParticleGroup)
    th = theta(pg)
    return pg.px .* cos.(th) .+ pg.py .* sin.(th)
end

function ptheta(pg::AbstractParticleGroup)
    th = theta(pg)
    return -pg.px .* sin.(th) .+ pg.py .* cos.(th)
end

Lz(pg::AbstractParticleGroup) = pg.x .* pg.py .- pg.y .* pg.px

gamma(pg::AbstractParticleGroup) = energy(pg) ./ mass(pg)
beta(pg::AbstractParticleGroup) = momentum(pg) ./ energy(pg)
beta_x(pg::AbstractParticleGroup) = pg.px ./ energy(pg)
beta_y(pg::AbstractParticleGroup) = pg.py ./ energy(pg)
beta_z(pg::AbstractParticleGroup) = pg.pz ./ energy(pg)

# These call normalized_particle_coordinate (defined in statistics.jl, resolved at call time)
x_bar(pg::AbstractParticleGroup) = normalized_particle_coordinate(pg, :x)
px_bar(pg::AbstractParticleGroup) = normalized_particle_coordinate(pg, :px)
y_bar(pg::AbstractParticleGroup) = normalized_particle_coordinate(pg, :y)
py_bar(pg::AbstractParticleGroup) = normalized_particle_coordinate(pg, :py)

Jx(pg::AbstractParticleGroup) = particle_amplitude(pg; plane=:x)
Jy(pg::AbstractParticleGroup) = particle_amplitude(pg; plane=:y)

# ---- set_charge! ----

function set_charge!(pg::AbstractParticleGroup, val)
    val > 0 || throw(ArgumentError("charge must be >0. This is used to weight the particles."))
    pg.weight .*= val / charge(pg)
end

# ---- Display ----

function Base.show(io::IO, pg::ParticleGroup)
    print(io, "ParticleGroup with $(length(pg)) particles (charge = $(charge(pg)) C)")
end

function Base.show(io::IO, pg::ParticleGroupView)
    print(io, "ParticleGroupView with $(length(pg)) particles (charge = $(charge(pg)) C)")
end

# ---- copy ----

"""
    copy(pg::AbstractParticleGroup)

Return a new `ParticleGroup` with copied arrays. Useful for preserving state
before in-place operations like `drift!`.
"""
function Base.copy(pg::AbstractParticleGroup)
    return ParticleGroup(
        copy(pg.x), copy(pg.px), copy(pg.y), copy(pg.py),
        copy(pg.z), copy(pg.pz), copy(pg.t), copy(pg.status),
        copy(pg.weight), pg.species, copy(pg.id), Dict{Symbol, AbstractVector}(k => copy(v) for (k, v) in pg.extra)
    )
end

# Resolve a stored column or evaluate a user-supplied column function.
_particle_values(pg, column::Symbol) = getproperty(pg, column)::AbstractVector

function _particle_values(pg, f)
    values = f(pg)::AbstractVector
    length(values) == length(pg) || throw(DimensionMismatch("Expected one value per particle"))
    return values
end
# ---- Coordinate checks ----

function in_z_coordinates(pg::AbstractParticleGroup)
    return all(==(pg.z[1]), pg.z)
end

function in_t_coordinates(pg::AbstractParticleGroup)
    return all(==(pg.t[1]), pg.t)
end

# ---- average_current ----

function average_current(pg::AbstractParticleGroup)
    dt = ptp(pg.t)
    if dt == 0
        dt = ptp(pg.z) / (mean(pg, beta_z) * C_LIGHT)
    end
    return charge(pg) / dt
end

# ---- Drift methods ----

function drift!(pg::AbstractParticleGroup, delta_t::Union{Number,AbstractVector})
    if delta_t isa AbstractVector && length(delta_t) != length(pg)
        throw(ArgumentError("delta_t length must match number of particles"))
    end
    E = energy(pg)
    pg.x .+= (pg.px ./ E) .* C_LIGHT .* delta_t
    pg.y .+= (pg.py ./ E) .* C_LIGHT .* delta_t
    pg.z .+= (pg.pz ./ E) .* C_LIGHT .* delta_t
    pg.t .+= delta_t
end

function drift_to_z!(pg::AbstractParticleGroup, z::Real)
    E = energy(pg)
    dt = (z .- pg.z) ./ ((pg.pz ./ E) .* C_LIGHT)
    drift!(pg, dt)
    pg.z .= z
end

function drift_to_z!(pg::AbstractParticleGroup)
    z = mean(pg, :z)
    drift_to_z!(pg, z)
end

function drift_to_t!(pg::AbstractParticleGroup, t::Real)
    dt = t .- pg.t
    drift!(pg, dt)
    pg.t .= t
end

function drift_to_t!(pg::AbstractParticleGroup)
    t = mean(pg, :t)
    drift_to_t!(pg, t)
end

# Indexing an owning group copies rows; indexing a view preserves the view.
Base.getindex(pg::AbstractParticleGroup, i::Integer) = pg[i:i]
Base.getindex(pg::ParticleGroup, rows::Union{Colon,AbstractVector{<:Integer}}) =
    particle_parts(pg, rows)
Base.getindex(pg::ParticleGroupView, rows::Union{Colon,AbstractVector{<:Integer}}) =
    view(pg, rows)

# ---- Arithmetic and comparison ----

function Base.:+(pg1::AbstractParticleGroup, pg2::AbstractParticleGroup)
    return join_particle_groups(pg1, pg2)
end

function Base.:(==)(pg1::AbstractParticleGroup, pg2::AbstractParticleGroup)
    if nameof(pg1.species) != nameof(pg2.species)
        return false
    end
    for key in (:x, :px, :y, :py, :z, :pz, :t, :status, :weight, :id)
        if getfield(pg1, key) != getfield(pg2, key)
            return false
        end
    end
    return pg1.extra == pg2.extra
end

function Base.isapprox(pg1::AbstractParticleGroup, pg2::AbstractParticleGroup; kwargs...)
    if nameof(pg1.species) != nameof(pg2.species)
        return false
    end
    for key in (:x, :px, :y, :py, :z, :pz, :t, :weight)
        if !isapprox(getfield(pg1, key), getfield(pg2, key); kwargs...)
            return false
        end
    end
    for key in (:status, :id)
        if getfield(pg1, key) != getfield(pg2, key)
            return false
        end
    end
    keys(pg1.extra) == keys(pg2.extra) || return false
    return all(isapprox(pg1.extra[k], pg2.extra[k]; kwargs...) for k in keys(pg1.extra))
end

# ---- Helper functions ----

"""
    slices(pg; by=nothing, n=40)

Sort the group by the selected column and return n roughly equal-count views.
When by is omitted, use z for a common-time group and t otherwise.
"""
function slices(pg::AbstractParticleGroup; by=nothing, n=40)
    column = isnothing(by) ? (in_t_coordinates(pg) ? :z : :t) : by
    order = sortperm(_particle_values(pg, column))
    count = length(order)
    1 <= n <= count || throw(ArgumentError("Cannot split $count particles into $n slices"))
    edges = round.(Int, range(1, count + 1, length=n + 1))
    return [view(pg, @view(order[edges[i]:edges[i+1]-1])) for i in 1:n]
end

function particle_parts(pg::ParticleGroup, idx)
    return ParticleGroup(
        pg.x[idx], pg.px[idx], pg.y[idx], pg.py[idx],
        pg.z[idx], pg.pz[idx], pg.t[idx], pg.status[idx],
        pg.weight[idx], pg.species, pg.id[idx], Dict{Symbol, AbstractVector}(k => v[idx] for (k, v) in pg.extra)
    )
end

function join_particle_groups(pgs::AbstractParticleGroup...)
    species = [pg.species for pg in pgs]
    species0 = species[1]
    @assert all(spe == species0 for spe in species) "species must be the same to join"

    x = reduce(vcat, (pg.x for pg in pgs))
    px = reduce(vcat, (pg.px for pg in pgs))
    y = reduce(vcat, (pg.y for pg in pgs))
    py = reduce(vcat, (pg.py for pg in pgs))
    z = reduce(vcat, (pg.z for pg in pgs))
    pz = reduce(vcat, (pg.pz for pg in pgs))
    t = reduce(vcat, (pg.t for pg in pgs))
    status = reduce(vcat, (pg.status for pg in pgs))
    weight = reduce(vcat, (pg.weight for pg in pgs))
    id = reduce(vcat, (pg.id for pg in pgs))
    all_extra = union((keys(pg.extra) for pg in pgs)...)
    common = intersect((keys(pg.extra) for pg in pgs)...)
    dropped = setdiff(all_extra, common)
    isempty(dropped) || @warn "Dropping optional columns not present in every group" dropped
    extra = Dict{Symbol, AbstractVector}(k => reduce(vcat, (pg.extra[k] for pg in pgs)) for k in common)

    return ParticleGroup(
        x, px, y, py, z, pz, t, status, weight, species0, id, extra
    )
end

# ---- Particle file I/O ----

"""
    ParticleGroup(h5::Union{HDF5.Group,HDF5.File}; include_time_offset=true)

Load a particle group from an open HDF5 file or group in openPMD-beamphysics format.
The group may point directly to a species subgroup or to a path containing
exactly one species. Positions are returned in m, momenta in eV/c, and time in s.
Set include_time_offset=false to read the time record without adding timeOffset.
"""
function ParticleGroup(h5::Union{HDF5.Group,HDF5.File}; include_time_offset=true)
    if haskey(attrs(h5), "basePath") && haskey(attrs(h5), "particlesPath")
        paths = particle_paths(h5)
        length(paths) == 1 || throw(ArgumentError("Expected one particle iteration, found $(length(paths))"))
        path = strip(only(paths), '/')
        if !isempty(path)
            h5 = h5[path]
        end
    end
    if !haskey(h5, "position")
        species = keys(h5)
        if length(species) != 1
            error("multiple species in particle paths: $(species)")
        end
        h5 = h5[species[1]]
    end

    attributes = attrs(h5)
    species = Species(_scalar_attr(attributes["speciesType"]))
    n_particle = Int(_scalar_attr(attributes["numParticles"]))
    total_charge = _scalar_attr(attributes["totalCharge"]) * _scalar_attr(attributes["chargeUnitSI"])

    x = particle_array(h5, "x")
    px = particle_array(h5, "px")
    y = particle_array(h5, "y")
    py = particle_array(h5, "py")
    z = particle_array(h5, "z")
    pz = particle_array(h5, "pz")
    t = particle_array(h5, "t"; include_offset=include_time_offset)

    if haskey(h5, "particleStatus")
        status = particle_array(h5, "particleStatus")
    else
        status = fill(1, n_particle)
    end

    if haskey(h5, "weight")
        weight = particle_array(h5, "weight")
        length(weight) == 1 && (weight = fill(only(weight), n_particle))
    else
        weight = fill(total_charge / n_particle, n_particle)
    end

    if haskey(h5, "id")
        id = particle_array(h5, "id")
    else
        id = collect(1:n_particle)
    end
    extra = Dict{Symbol, AbstractVector}()
    for key in OPTIONAL_PARTICLE_COLUMNS
        alias = string(key)
        if haskey(h5, COMPONENT_FROM_ALIAS[alias])
            extra[key] = particle_array(h5, alias)
        end
    end
    # A file can mix floating types across its records; a group uses one type.
    T = promote_type(map(eltype, (x, px, y, py, z, pz, t, weight))...)
    x, px, y, py, z, pz, t, weight = map(v -> convert(Vector{T}, v), (x, px, y, py, z, pz, t, weight))
    return ParticleGroup(; x, px, y, py, z, pz, t, status, weight, species, id, extra)
end

"""
    ParticleGroup(file::AbstractString; include_time_offset=true)

Load particles into a ParticleGroup from an HDF5 file in openPMD format.
"""
function ParticleGroup(file::AbstractString; include_time_offset=true)
    h5open(file, "r") do h5
        return ParticleGroup(h5; include_time_offset)
    end
end
