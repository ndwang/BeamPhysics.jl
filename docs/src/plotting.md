# Plotting

Plotting functionality is available as a package extension when
[Plots.jl](https://github.com/JuliaPlots/Plots.jl) is loaded.

## Setup

```julia
using BeamPhysics
using Plots
```

## Particle Plots

```julia
# 1D density histogram of a particle property
density_plot(pg, :x)
density_plot(pg, :z; bins=100)

# Marginal plot: 2D density with 1D projections on margins
marginal_plot(pg, :x, :y)
marginal_plot(pg, :z, energy)

# Slice plot: statistics along slices
slice_plot(pg, norm_emit_x, norm_emit_y)
slice_plot(pg, norm_emit_x; slice_key=:z)

# Combined density + slice
density_and_slice_plot(pg, :z, :x)
```

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
