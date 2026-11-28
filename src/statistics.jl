_masked_group(pg, ::Nothing) = pg

function _masked_group(pg, mask::AbstractVector{Bool})
    length(mask) == length(pg) || throw(DimensionMismatch("mask must have one entry per particle"))
    rows = findall(mask)
    isempty(rows) && throw(ArgumentError("mask selects no particles"))
    return view(pg, rows)
end

"""
    mean(pg, column; mask=nothing)
    std(pg, column; mask=nothing)

Weighted mean and population standard deviation. Pass a symbol naming a stored
column or a function returning one value per particle. A Boolean mask selects
particles before evaluating the column. By default, all particles contribute.
"""
function mean(pg::AbstractParticleGroup, column; mask=nothing)
    pg = _masked_group(pg, mask)
    values = _particle_values(pg, column)
    length(values) == 1 && return only(values)
    return mean(values, weights(pg.weight))
end

function std(pg::AbstractParticleGroup, column; mask=nothing)
    pg = _masked_group(pg, mask)
    values = _particle_values(pg, column)
    length(values) == 1 && return zero(eltype(values))
    w = weights(pg.weight)
    center = mean(values, w)
    return sqrt(mean((values .- center).^2, w))
end

Base.minimum(pg::AbstractParticleGroup, column) = minimum(_particle_values(pg, column))
Base.maximum(pg::AbstractParticleGroup, column) = maximum(_particle_values(pg, column))

"""
    cov(pg, first, second; mask=nothing)

Weighted population covariance between two particle columns. Each selector is a
stored column symbol or a function returning one value per particle. Pass a
Boolean mask to compute the covariance over selected particles.
"""
function cov(pg::AbstractParticleGroup, first, second; mask=nothing)
    pg = _masked_group(pg, mask)
    x = _particle_values(pg, first)
    y = _particle_values(pg, second)
    w = weights(pg.weight)
    x0, y0 = mean(x, w), mean(y, w)
    return mean((x .- x0) .* (y .- y0), w)
end

"""
    covariance_matrix(pg, columns...; mask=nothing)

Weighted population covariance matrix, with rows and columns in selector order.
Selectors and the optional Boolean mask are the same as for `cov`.
"""
function covariance_matrix(pg::AbstractParticleGroup, columns...; mask=nothing)
    pg = _masked_group(pg, mask)
    data = hcat(map(column -> _particle_values(pg, column), columns)...)
    return cov(data, weights(pg.weight), 1; corrected=false)
end

function delta(pg::AbstractParticleGroup, column)
    values = _particle_values(pg, column)
    center = length(values) == 1 ? only(values) : mean(values, weights(pg.weight))
    return values .- center
end

ptp(x::AbstractArray) = maximum(x) - minimum(x)
ptp(pg::AbstractParticleGroup, column) = ptp(_particle_values(pg, column))

function _momentum_column(plane::Symbol)
    plane in (:x, :y, :z) || throw(ArgumentError("Unknown plane: $plane"))
    return Symbol(:p, plane)
end

"""
    norm_emit_calc(pg, planes=(:x,))

Compute normalized emittance in the selected planes using mechanical momenta.
For N planes, the units are m^N.
"""
function norm_emit_calc(pg::AbstractParticleGroup, planes=(:x,))
    columns = map(plane -> (plane, _momentum_column(plane)), planes)
    S = covariance_matrix(pg, Iterators.flatten(columns)...)
    return sqrt(det(S)) / mass(pg)^length(planes)
end

norm_emit_x(pg::AbstractParticleGroup) = norm_emit_calc(pg, (:x,))
norm_emit_y(pg::AbstractParticleGroup) = norm_emit_calc(pg, (:y,))
norm_emit_4d(pg::AbstractParticleGroup) = norm_emit_calc(pg, (:x, :y))

"""
    twiss_calc(sigma_mat2)

Calculate Twiss parameters from the 2D sigma matrix (covariance matrix):
```
sigma_mat = <x,x>   <x, p>
            <p, x>  <p, p>
```

This is a simple calculation. Makes no assumptions about units.

Returns a NamedTuple with fields: alpha, beta, gamma, emit
"""
function twiss_calc(sigma_mat2::AbstractMatrix{T}) where T
    @assert size(sigma_mat2) == (2, 2) "Bad shape: $(size(sigma_mat2)). This should be (2,2)"

    emit = sqrt(det(sigma_mat2))
    return (;
        alpha = -sigma_mat2[1, 2] / emit,
        beta  =  sigma_mat2[1, 1] / emit,
        gamma =  sigma_mat2[2, 2] / emit,
        emit  =  emit
    )
end

"""
    twiss_dispersion(sigma3)

Twiss and Dispersion calculation from a 3x3 sigma (covariance) matrix from particles
x, p, delta

Returns a NamedTuple with fields: alpha, beta, gamma, emit, eta, etap
"""
function twiss_dispersion(sigma3::Matrix)
    delta2 = sigma3[3, 3]
    xd = sigma3[1, 3]
    pd = sigma3[2, 3]

    eb = sigma3[1, 1] - xd^2 / delta2
    eg = sigma3[2, 2] - pd^2 / delta2
    ea = -sigma3[1, 2] + xd * pd / delta2

    emit = sqrt(eb * eg - ea^2)

    return (;
        alpha = ea / emit,
        beta  = eb / emit,
        gamma = eg / emit,
        emit  = emit,
        eta   = xd / delta2,
        etap  = pd / delta2
    )
end

"""
    twiss(pg; plane=:x, fraction=1, p0c=nothing)

Calculate dispersion-corrected Twiss parameters in plane :x or :y. Return a
named tuple with alpha, beta, gamma, emit, eta, etap, and norm_emit.

If p0c is omitted, use the weighted mean momentum of the selected particles.
For fraction < 1, select the lowest-action particles by count, not by charge;
the selected particles still contribute with their weights.
"""
function twiss(pg::AbstractParticleGroup; plane=:x, fraction=1, p0c=nothing)
    plane in (:x, :y) || throw(ArgumentError("Twiss plane must be :x or :y"))
    if fraction < 1
        action = particle_amplitude(pg; plane)
        order = sortperm(action)
        pg = view(pg, @view(order[1:trunc(Int, fraction * length(pg))]))
    end

    p = momentum(pg)
    if isnothing(p0c)
        p0c = length(p) == 1 ? only(p) : mean(p, weights(pg.weight))
    end

    x = getproperty(pg, plane)
    p_transverse = getproperty(pg, _momentum_column(plane)) ./ p0c
    sigma3 = cov(hcat(x, p_transverse, p ./ p0c), weights(pg.weight); corrected=false)
    tw = twiss_dispersion(sigma3)
    return (; tw..., norm_emit = tw.emit * p0c / mass(pg))
end

"""
    A_mat_calc(beta, alpha)

Returns the 1D normal form matrix from twiss parameters beta and alpha

    A =   sqrt(beta)         0
         -alpha/sqrt(beta)   1/sqrt(beta)
"""
function A_mat_calc(beta, alpha)
    a11 = sqrt(beta)
    a22 = 1 / a11
    a21 = -alpha / a11
    return [a11 0; a21 a22]
end

"""
    A_inverse_mat_calc(beta, alpha)

Returns the inverse of the 1D normal form matrix from twiss parameters beta and alpha

    A^-1 =  1/sqrt(beta)     0
            alpha/sqrt(beta) sqrt(beta)
"""
function A_inverse_mat_calc(beta, alpha)
    a11 = sqrt(beta)
    a22 = 1 / a11
    a21 = -alpha / a11
    return [a22 0; -a21 a11]
end

"""
    amplitude_calc(x, p; beta=1, alpha=0)

Simple amplitude calculation of position and momentum coordinates
relative to twiss beta and alpha.

J = (gamma x^2 + 2 alpha x p + beta p^2)/2
"""
function amplitude_calc(x, p; beta=1, alpha=0)
    return (1 + alpha^2) / beta / 2 .* x.^2 .+ alpha .* x .* p .+ beta / 2 .* p.^2
end

function _coordinates_and_twiss(pg, plane; mass_normalize, twiss)
    x = getproperty(pg, plane)
    p = getproperty(pg, _momentum_column(plane))
    if mass_normalize
        p = p ./ mass(pg)
    end
    if isnothing(twiss)
        sigma_mat2 = cov(hcat(x, p), weights(pg.weight), 1; corrected=false)
        twiss = twiss_calc(sigma_mat2)
    end
    return x, p, twiss
end

"""
    particle_amplitude(particle_group; plane=:x, twiss=nothing, mass_normalize=true)

Returns the action J = (x_bar^2 + p_bar^2)/2 for each particle in the given plane.
"""
function particle_amplitude(particle_group::AbstractParticleGroup; plane=:x, twiss=nothing, mass_normalize=true)
    x, p, twiss = _coordinates_and_twiss(particle_group, plane; mass_normalize, twiss)
    return amplitude_calc(x, p; beta=twiss.beta, alpha=twiss.alpha)
end

"""
    normalized_particle_coordinate(particle_group, key::Symbol; twiss=nothing, mass_normalize=true)

Returns a single normalized coordinate array from a ParticleGroup.
"""
function normalized_particle_coordinate(particle_group::AbstractParticleGroup, key::Symbol; twiss=nothing, mass_normalize=true)
    plane = get((px=:x, py=:y, pz=:z), key, key)
    x, p, twiss = _coordinates_and_twiss(particle_group, plane; mass_normalize, twiss)
    A_inv = A_inverse_mat_calc(twiss.beta, twiss.alpha)

    if key != plane
        return A_inv[2, 1] .* x .+ A_inv[2, 2] .* p
    else
        return A_inv[1, 1] .* x
    end
end
