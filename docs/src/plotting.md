# Plotting

Plotting functionality is available as a package extension when
[Plots.jl](https://github.com/JuliaPlots/Plots.jl) is loaded.

## Setup

```julia
using BeamPhysics
using Plots
```

## Particle Plots

The Quick Start bunch has a common z, so its longitudinal distribution is plotted against time.

```julia
# 1D density histogram of a particle property
density_plot(pg, :x)
density_plot(pg, :t; bins=100)

# Marginal plot: 2D density with 1D projections on margins
marginal_plot(pg, :x, :y)
marginal_plot(pg, :t, energy)

# Slice plot: statistics along slices
slice_plot(pg, norm_emit_x, norm_emit_y)
slice_plot(pg, norm_emit_x; slice_key=:t)

# A custom function returns one value per slice
width_mm(g) = 1000 * std(g, :x)
slice_plot(pg, "Beam width (mm)" => width_mm)

# Time–energy heatmap with slice statistics on the right axis
density_and_slice_plot(pg, :t, energy, norm_emit_x, norm_emit_y)
density_and_slice_plot(pg, :t, energy, "Beam width (mm)" => width_mm)
```

In the combined plot, the first coordinate also selects the slicing direction.
Both plots accept one or more slice statistics as positional arguments. If either
plot encounters an unrecognized name or label, it displays all slice statistics in
their supplied units. Include those units in custom labels.

See [`density_plot`](@ref), [`marginal_plot`](@ref), [`slice_plot`](@ref),
[`density_and_slice_plot`](@ref).

## Field Mesh Plots

Both plotting functions use the physical fields at the mesh's current scale
and RF phase. After changing a setting or editing the profile, call the plot
function again to see the result. See [Scaling an RF field map](@ref).

```julia
# 2D field plot (slice through the mesh)
fieldmesh_plot(fm)                             # auto-detects :B or :E
fieldmesh_plot(fm; component=:Bz)              # specific component
fieldmesh_plot(fm; component=:Bz, mirror=true) # mirror for cylindrical

# On-axis field profile
fieldmesh_plot_onaxis(fm)
```

On a rectangular mesh, `fieldmesh_plot(fm; component=:Bz, z=0.4)` shows
x horizontally and y vertically. The field is linearly interpolated at
the requested z coordinate in metres. Use `x=value` or `y=value` for a
longitudinal slice. The default slice is at y = 0.

For an RF field, `fieldmesh_plot(fm; component=:Ez, transform=imag)` displays
the imaginary part. Use `transform=abs` for the component amplitude. The default
is `real`. The transform is applied at the grid points before interpolating the slice.
On-axis RF profiles label the real and imaginary traces separately.
On rectangular meshes these profiles interpolate the field at x = y = 0,
including when the axis lies between grid points.

See [`fieldmesh_plot`](@ref), [`fieldmesh_plot_onaxis`](@ref).
