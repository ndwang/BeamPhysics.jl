module BeamPhysicsPlotsExt

using Plots
using BeamPhysics
using BeamPhysics: pg_units, plottable_array, nice_array, nice_scale_prefix,
    mathlabel, ptp, limits, in_t_coordinates, slices, covariance_matrix, momentum, _particle_values

import BeamPhysics: slice_plot, density_plot, marginal_plot, density_and_slice_plot
import BeamPhysics: fieldmesh_plot, fieldmesh_plot_onaxis

using StatsBase: weights, fit, Histogram
import Statistics

# Default colormaps
const CMAP0 = :viridis
const CMAP1 = :plasma

# Sigma ellipse for 2x2 covariance matrix (1-sigma)
function _sigma_ellipse(sigma_mat2; n=100)
    a, b = sigma_mat2[1, 1], sigma_mat2[1, 2]
    d = sigma_mat2[2, 2]
    tr = (a + d) / 2
    disc = sqrt((a - d)^2 / 4 + b^2)
    λ1, λ2 = tr + disc, tr - disc
    angle = atan(2b, a - d) / 2
    t = range(0, 2π, length=n)
    ca, sa = cos(angle), sin(angle)
    r1, r2 = sqrt(max(λ1, 0)), sqrt(max(λ2, 0))
    x = r1 .* ca .* cos.(t) .- r2 .* sa .* sin.(t)
    y = r1 .* sa .* cos.(t) .+ r2 .* ca .* sin.(t)
    return x, y
end

_plot_key(column::Symbol) = string(column)
_plot_key(f::Function) = string(nameof(f))
_plot_key(::typeof(momentum)) = "p"
_plot_key(diagnostic::Pair) = string(first(diagnostic))
_plot_function(f) = f
_plot_function(diagnostic::Pair) = last(diagnostic)

function _slice_plot_data(pg, by, n, diagnostics)
    groups = slices(pg; by, n)
    center = [Statistics.mean(s, by) for s in groups]
    density = [charge(s) / ptp(s, by) for s in groups]
    values = map(diagnostic -> map(_plot_function(diagnostic), groups), diagnostics)
    return center, density, values
end

"""
    slice_plot(particle_group, diagnostics...; n_slice=40, slice_key=nothing, xlim=nothing, ylim=nothing, tex=true, kwargs...)

Complete slice plotting routine. Will plot the density of the slice key on the right axis.

# Arguments
- `particle_group`: ParticleGroup object to plot
- Diagnostics: pass functions such as norm_emit_x, or a label => function pair.
- `n_slice`: Number of slices (default: 40)
- `slice_key`: Dimension to slice in (default: automatically determined)
- `xlim`: Manual setting of x-axis limits
- `ylim`: Manual setting of y-axis limits
- `tex`: Use TEX for labels (default: true)

# Returns
- `p`: Plots.Plot object
"""
function slice_plot(particle_group, diagnostics...; n_slice=40, slice_key=nothing, xlim=nothing, ylim=nothing, tex=true, kwargs...)
    by = isnothing(slice_key) ? (in_t_coordinates(particle_group) ? :z : :t) : slice_key
    keys = map(_plot_key, diagnostics)
    x, density, values = _slice_plot_data(particle_group, by, n_slice, diagnostics)
    slice_key = _plot_key(by)
    # X-axis
    x, f1, p1, xmin, xmax = plottable_array(x; lim=xlim)
    ux = p1 * string(pg_units(slice_key))

    # Y-axis
    ulist = [pg_units(k) for k in keys]
    uy = ulist[1]
    if !all(u == uy for u in ulist)
        error("Incompatible units: $ulist")
    end
    uy_str = string(uy)

    ymin = minimum(minimum, values)
    ymax = maximum(maximum, values)

    _, f2, p2, ymin, ymax = plottable_array([ymin, ymax]; lim=ylim)
    uy_str = p2 * uy_str

    # Calculate plot limits
    plot_xlim = isnothing(xlim) ? (xmin/f1, xmax/f1) : xlim
    plot_ylim = isnothing(ylim) ? (ymin/f2, ymax/f2) : ylim

    # Create plot
    p = plot(; xlim=plot_xlim, ylim=plot_ylim, kwargs...)

    # Main curves
    color = length(keys) == 1 ? :black : nothing

    for (k, y) in zip(keys, values)
        label = mathlabel(k, units=uy_str, tex=tex)
        plot!(p, x, y ./ f2, label=label, color=color)
    end

    if length(keys) > 1
        plot!(p, legend=true)
    end

    # Density on right axis
    y2, _, prey2, _, _ = plottable_array(density; lim=nothing)

    # Convert to Amps if possible
    y2_units = "C/" * string(pg_units(slice_key))
    if y2_units == "C/s"
        y2_units = "A"
    end
    y2_units = prey2 * y2_units

    # Labels
    labelx = mathlabel(slice_key, units=ux, tex=tex)
    labely = mathlabel(keys..., units=uy, tex=tex)
    labely2 = mathlabel("density", units=y2_units, tex=tex)

    xlabel!(p, labelx)
    ylabel!(p, labely)

    # Right axis plot
    p2 = twinx(p)
    ylabel!(p2, labely2)
    plot!(p2, x, y2, fillrange=0, fillalpha=0.2, color=:gray, w=0)

    return p
end

"""
    density_plot(particle_group, key=:x; bins=nothing, xlim=nothing, tex=true, kwargs...)

1D density plot.

# Arguments
- `particle_group`: ParticleGroup object to plot
- `key`: Key to plot (default: :x)
- `bins`: Number of bins (default: automatically determined)
- `xlim`: Manual setting of x-axis limits
- `tex`: Use TEX for labels (default: true)

# Returns
- `p`: Plots.Plot object
"""
function density_plot(particle_group, key=:x; bins=nothing, xlim=nothing, tex=true, kwargs...)
    if isnothing(bins)
        n = length(particle_group)
        bins = max(1, floor(Int, sqrt(n/4)))
    end

    # Scale to nice units
    x, f1, p1, xmin, xmax = plottable_array(_particle_values(particle_group, key); lim=xlim)
    w = particle_group.weight
    u1 = string(pg_units(_plot_key(key)))
    ux = p1 * u1

    # Create label
    labelx = mathlabel(_plot_key(key), units=ux, tex=tex)

    # Create histogram
    hist = fit(Histogram, x, weights(w), nbins=bins)
    hist_x = hist.edges[1][1:end-1] .+ diff(hist.edges[1])./2
    hist_width = diff(hist.edges[1])
    hist_y, hist_f, hist_prefix = nice_array(hist.weights ./ hist_width)

    plot_xlim = isnothing(xlim) ? (xmin/f1, xmax/f1) : xlim
    # Create plot
    p = bar(hist_x, hist_y, bar_width=hist_width, lw=0,xlim=plot_xlim, color=:gray; kwargs...)

    # Set labels
    xlabel!(p, labelx)
    if u1 == "s"
        _, hist_prefix = nice_scale_prefix(hist_f/f1)
        ylabel!(p, "$(hist_prefix)A")
    else
        ylabel!(p, "$(hist_prefix)C/$ux")
    end

    return p
end

"""
    marginal_plot(particle_group, key1=:t, key2=momentum; bins=nothing, xlim=nothing, ylim=nothing, tex=true, ellipse=false, kwargs...)

Density plot and projections.

# Arguments
- `particle_group`: ParticleGroup object to plot
- `key1`: Key for x-axis (default: :t)
- `key2`: Key for y-axis (default: momentum)
- `bins`: Number of bins (default: automatically determined)
- `xlim`: Manual setting of x-axis limits
- `ylim`: Manual setting of y-axis limits
- `tex`: Use TEX for labels (default: true)
- `ellipse`: Plot sigma matrix ellipse (default: false)

# Returns
- `p`: Plots.Plot object
"""
function marginal_plot(particle_group, key1=:t, key2=momentum; bins=nothing, xlim=nothing, ylim=nothing, tex=true, ellipse=false, kwargs...)
    if isnothing(bins)
        n = length(particle_group)
        bins = max(1, floor(Int, sqrt(n/4)))
    end

    # Get data
    x = _particle_values(particle_group, key1)
    y = _particle_values(particle_group, key2)

    if length(x) == 1
        bins = 100

        if isnothing(xlim)
            x0 = x[1]
            if isapprox(x0, 0.0)
                xlim = (-1, 1)
            else
                xlim = sort([0.9*x0, 1.1*x0])
            end
        end
        if isnothing(ylim)
            y0 = y[1]
            if isapprox(y0, 0.0)
                ylim = (-1, 1)
            else
                ylim = sort([0.9*y0, 1.1*y0])
            end
        end
    end

    # Scale to nice units
    x, f1, p1, xmin, xmax = plottable_array(x; lim=xlim)
    y, f2, p2, ymin, ymax = plottable_array(y; lim=ylim)
    w = particle_group.weight

    u1 = string(pg_units(_plot_key(key1)))
    u2 = string(pg_units(_plot_key(key2)))
    ux = p1 * u1
    uy = p2 * u2

    # Create labels
    labelx = mathlabel(_plot_key(key1), units=ux, tex=tex)
    labely = mathlabel(_plot_key(key2), units=uy, tex=tex)

    # Calculate plot limits
    plot_xlim = isnothing(xlim) ? (xmin/f1, xmax/f1) : xlim
    plot_ylim = isnothing(ylim) ? (ymin/f2, ymax/f2) : ylim

    # Create plot
    p = plot(; xlim=plot_xlim, ylim=plot_ylim, kwargs...)

    if all(isnan, x)
        annotate!(p, 0.5, 0.5, text("$key1 is all NaN", :center))
        return p
    end
    if all(isnan, y)
        annotate!(p, 0.5, 0.5, text("$key2 is all NaN", :center))
        return p
    end

    # Main plot
    if length(x) == 1
        scatter!(p, x, y, xlim=plot_xlim, ylim=plot_ylim)
    else
        histogram2d!(p, x, y, weights(w), bins=bins, color=CMAP0, clims=(1e-20, Inf), xlim=plot_xlim, ylim=plot_ylim, colorbar=false)
    end

    # Add ellipse if requested
    if ellipse
        sigma_mat2 = covariance_matrix(particle_group, key1, key2)
        x_ellipse, y_ellipse = _sigma_ellipse(sigma_mat2)
        x_ellipse .+= Statistics.mean(particle_group, key1)
        y_ellipse .+= Statistics.mean(particle_group, key2)
        plot!(p, x_ellipse ./ f1, y_ellipse ./ f2, color=:red, label="")
    end

    # Add histograms
    hist_x = fit(Histogram, x, weights(w), nbins=bins)
    hist_y = fit(Histogram, y, weights(w), nbins=bins)

    # Top histogram
    hist_x_centers = hist_x.edges[1][1:end-1] .+ diff(hist_x.edges[1])./2
    hist_x_width = diff(hist_x.edges[1])
    hist_x_y, hist_x_f, hist_x_prefix = nice_array(hist_x.weights ./ hist_x_width)

    p_top = bar(hist_x_centers, hist_x_y, bar_width=hist_x_width, lw=0,color=:grey, xlim=plot_xlim,
                xticks=false, xlabel="", legend=false)
    if u1 == "s"
        _, hist_x_prefix = nice_scale_prefix(hist_x_f/f1)
        ylabel!(p_top, "$(hist_x_prefix)A")
    else
        ylabel!(p_top, "$(hist_x_prefix)" * mathlabel("C/$ux"))
    end

    # Side histogram
    hist_y_centers = hist_y.edges[1][1:end-1] .+ diff(hist_y.edges[1])./2
    hist_y_width = diff(hist_y.edges[1])
    hist_y_y, hist_y_f, hist_y_prefix = nice_array(hist_y.weights ./ hist_y_width)

    p_side = bar(hist_y_centers, hist_y_y, bar_width=hist_y_width, lw=0,color=:grey, ylim=plot_ylim, orientation=:h,
                 yticks=false, ylabel="", legend=false)
    xlabel!(p_side, "$(hist_y_prefix)" * mathlabel("C/$uy"))

    # Set labels and limits
    xlabel!(p, labelx)
    ylabel!(p, labely)

    # Combine plots
    layout = @layout [a{0.3h} _; b c{0.3w}]
    return plot(p_top, p, p_side, layout=layout)
end

"""
    density_and_slice_plot(particle_group, key1=:t, key2=momentum; diagnostics=(norm_emit_x, norm_emit_y), bins=100, n_slice=30, tex=true)

Combined density and slice plot.

# Arguments
- `particle_group`: ParticleGroup object to plot
- `key1`: Key for x-axis (default: :t)
- `key2`: Key for y-axis (default: momentum)
- Diagnostics: functions to plot (default: normalized x and y emittance).
- `bins`: Number of bins (default: 100)
- `n_slice`: Number of slices (default: 30)
- `tex`: Use TEX for labels (default: true)

# Returns
- `p`: Plots.Plot object
"""
function density_and_slice_plot(particle_group, key1=:t, key2=momentum; diagnostics=(norm_emit_x, norm_emit_y), bins=100, n_slice=30, tex=true)
    # Scale to nice units
    x, f1, p1, xmin, xmax = plottable_array(_particle_values(particle_group, key1))
    y, f2, p2, ymin, ymax = plottable_array(_particle_values(particle_group, key2))
    w = particle_group.weight

    u1 = string(pg_units(_plot_key(key1)))
    u2 = string(pg_units(_plot_key(key2)))
    ux = p1 * u1
    uy = p2 * u2

    # Create labels
    labelx = mathlabel(_plot_key(key1), units=ux, tex=tex)
    labely = mathlabel(_plot_key(key2), units=uy, tex=tex)

    # Create plot
    p = plot()

    xlabel!(p, labelx)
    ylabel!(p, labely)

    # Create 2D histogram
    H = fit(Histogram, (x, y), weights(w), nbins=bins)
    heatmap!(p, H.edges[1], H.edges[2], H.weights', color=CMAP0, clims=(1e-16, Inf))

    centers, density, values = _slice_plot_data(particle_group, key1, n_slice, diagnostics)

    # Add slice statistics
    p2 = twinx(p)
    x2 = centers ./ f1
    stat_keys = map(_plot_key, diagnostics)
    ulist = map(pg_units, stat_keys)

    max2 = maximum(ptp, values)
    f3, p3 = nice_scale_prefix(max2)

    u2 = ulist[1]
    @assert all(u == u2 for u in ulist)
    u2 = p3 * string(u2)
    labely2 = mathlabel(stat_keys..., units=u2, tex=tex)

    for (k, y) in zip(stat_keys, values)
        label = mathlabel(k, units=u2, tex=tex)
        plot!(p2, x2, y ./ f3, label=label)
    end

    plot!(p2, legend=true)
    ylabel!(p2, labely2)

    # Add density
    y2 = density .* max2 ./ maximum(density) ./ f3 ./ 2
    plot!(p2, x2, y2, fillrange=0, fillalpha=0.2, color=:gray, w=0)
    return p
end


# ── FieldMesh plotting ─────────────────────────────────────────────────────

function _field_label(component, prefix, transform=nothing)
    name = string(component)
    label = mathlabel(name; units=prefix * string(pg_units(name)))
    transform === real && return "Re " * label
    transform === imag && return "Im " * label
    transform === abs && return "abs " * label
    return label
end

"""
    fieldmesh_plot(fm::FieldMesh; component, transform=real, mirror=false, cmap=:plasma, kwargs...)

2D heatmap of a FieldMesh field component.

For cylindrical geometry, plots the (r, z) plane at the first theta sample. Z is on the x-axis,
r on the y-axis (convention for accelerators).

For rectangular geometry, plots the (x, z) plane at y=0 by default. Pass
`x=value`, `y=value`, or `z=value` to interpolate at a coordinate within the grid.

# Arguments
- `component`: A native component such as `:Bz` or `:Ez`, or the vector magnitude
  `:B` or `:E`. Defaults to `:B` for a purely magnetic map and `:E` otherwise.
- `transform`: Apply `real` (default), `imag`, or `abs` to the field samples
  before interpolating a slice.
- `mirror`: For cylindrical only. If `true`, mirrors data to show negative-r side.
  `Br` and `Er` components are sign-flipped on the mirrored side.
- `cmap`: Colormap (default `:plasma`)
- `kwargs...`: Forwarded to `Plots.heatmap`

# Returns
`Plots.Plot`
"""
function fieldmesh_plot(fm::FieldMesh;
    component::Symbol = is_pure_magnetic(fm) ? :B : :E,
    transform = real,
    mirror::Bool = false,
    cmap = :plasma,
    kwargs...
)
    transform in (real, imag, abs) || throw(ArgumentError("transform must be real, imag, or abs"))
    values = fm[String(component)]

    if fm.geometry == "cylindrical"
        return _fieldmesh_plot_cylindrical(fm, component, values, transform; mirror, cmap, kwargs...)
    elseif fm.geometry == "rectangular"
        return _fieldmesh_plot_rectangular(fm, component, values, transform; cmap, kwargs...)
    else
        throw(ArgumentError("Unsupported geometry: $(fm.geometry)"))
    end
end

function _fieldmesh_plot_cylindrical(fm::FieldMesh, component, values, transform;
    mirror::Bool = false, cmap = :plasma, kwargs...)

    r = coord_vec(fm, "r")
    z = coord_vec(fm, "z")

    dat = transform.(@view values[:, 1, :])

    if mirror
        if r[1] != 0.0
            throw(ArgumentError("mirror=true only available when r=0 is at grid origin"))
        end
        # Sign convention: radial components flip sign, axial/magnitude do not
        sign = component in (:Br, :Er) && transform !== abs ? -1 : 1
        dat = vcat(sign .* reverse(dat; dims=1), dat)
        r = vcat(-reverse(r), r)
    end

    dat_scaled, fac, prefix = nice_array(dat)
    clabel = _field_label(component, prefix, transform === real && eltype(values) <: Real ? nothing : transform)
    xlabel = "\$z\$ (m)"
    ylabel = "\$r\$ (m)"

    return heatmap(z, r, dat_scaled;
        xlabel=xlabel, ylabel=ylabel, colorbar_title=clabel,
        color=cmap, kwargs...)
end

function _fieldmesh_plot_rectangular(fm::FieldMesh, component, values, transform;
    cmap = :plasma, kwargs...)

    # Pop slice coordinate from kwargs; default to y=0
    kw = Dict(kwargs)
    slice_axis, slice_val = :y, 0.0
    for ax in Symbol.(fm.axis_labels)
        if haskey(kw, ax)
            slice_axis = ax
            slice_val = Float64(pop!(kw, ax))
            break
        end
    end

    x, y, z = coord_vecs(fm)

    field = BeamPhysics.interpolator(fm, transform.(values))

    if slice_axis == :y
        dat = [field(xi, slice_val, zi) for xi in x, zi in z]
        hax, vax = z, x
        xlabel, ylabel = "\$z\$ (m)", "\$x\$ (m)"
    elseif slice_axis == :x
        dat = [field(slice_val, yi, zi) for yi in y, zi in z]
        hax, vax = z, y
        xlabel, ylabel = "\$z\$ (m)", "\$y\$ (m)"
    else  # z
        dat = [field(xi, yi, slice_val) for yi in y, xi in x]
        hax, vax = x, y
        xlabel, ylabel = "\$x\$ (m)", "\$y\$ (m)"
    end

    dat_scaled, fac, prefix = nice_array(dat)
    clabel = _field_label(component, prefix, transform === real && eltype(values) <: Real ? nothing : transform)

    return heatmap(hax, vax, dat_scaled;
        xlabel=xlabel, ylabel=ylabel, colorbar_title=clabel,
        color=cmap, kw...)
end

"""
    fieldmesh_plot_onaxis(fm::FieldMesh; kwargs...)

1D on-axis field plot for a FieldMesh.

For cylindrical geometry, plots `Bz` and/or `Ez` at r=0.
For rectangular geometry, plots all present field components along z at x=y=0.
RF profiles show their nonzero real and imaginary quadratures.

# Arguments
- `kwargs...`: Forwarded to `Plots.plot`

# Returns
`Plots.Plot`
"""
function fieldmesh_plot_onaxis(fm::FieldMesh; kwargs...)
    cylindrical = fm.geometry == "cylindrical"
    components = cylindrical ? (:Ez, :Bz) : (:Ex, :Ey, :Ez, :Bx, :By, :Bz)
    # The cylindrical axis uses r=0 at the map's first theta sample.
    second_coord = cylindrical ? fm.grid_origin[2] : 0.0
    z = coord_vec(fm, "z")
    p = plot(; xlabel="\$z\$ (m)", kwargs...)
    plotted_electric = false
    for component in components
        haskey(fm.components, BeamPhysics.COMPONENT_FROM_ALIAS[String(component)]) || continue
        field = BeamPhysics.interpolator(fm, String(component))
        values = [field(0, second_coord, zi) for zi in z]
        all(iszero, values) && continue
        if cylindrical
            axis = component === :Bz && plotted_electric ? twinx(p) : p
            color = component === :Bz ? :blue : :black
            _plot_field_profile!(axis, z, values, component; rf=!is_static(fm), color, label_axis=true)
            plotted_electric = component === :Ez
        else
            _plot_field_profile!(p, z, values, component; rf=!is_static(fm))
        end
    end
    return p
end

function _plot_field_profile!(p, z, values, component; rf, color=:auto, label_axis=false)
    re = real.(values)
    if rf
        im_part = imag.(values)
        factor, prefix = nice_scale_prefix(maximum(abs, values))
        if any(!iszero, re) || all(iszero, im_part)
            plot!(p, z, re ./ factor; label=_field_label(component, prefix, real), color)
        end
        if any(!iszero, im_part)
            plot!(p, z, im_part ./ factor; label=_field_label(component, prefix, imag),
                  color, linestyle=:dash)
        end
    else
        scaled, _, prefix = nice_array(re)
        plot!(p, z, scaled; label=_field_label(component, prefix), color)
    end
    label_axis && ylabel!(p, _field_label(component, prefix))
    return p
end

end # module
