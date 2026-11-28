# Field Mesh

## Overview

`FieldMesh` stores gridded electromagnetic field data in openPMD-beamphysics format.
It supports both rectangular (`x, y, z`) and cylindrical (`r, theta, z`) geometries
and both static and RF fields.

See [`FieldMesh`](@ref) in the [API Reference](@ref api) for full details.

## Loading and Saving

```julia
using BeamPhysics

# Load from file (must contain exactly one field path)
fm = FieldMesh("solenoid.h5")

# Save to file
write_fieldmesh("output.h5", fm)
```

See [`FieldMesh`](@ref), [`write_fieldmesh`](@ref).

## Accessing Fields and Coordinates

Field components are physical arrays in the mesh's native basis. The file's
scale and RF phase are already applied. Each array has shape `fm.grid_size`
and dimension order `fm.axis_labels`.

| Mesh geometry | Electric field [V/m] | Magnetic field [T] |
|---|---|---|
| Rectangular | `fm.Ex`, `fm.Ey`, `fm.Ez` | `fm.Bx`, `fm.By`, `fm.Bz` |
| Cylindrical | `fm.Er`, `fm.Etheta`, `fm.Ez` | `fm.Br`, `fm.Btheta`, `fm.Bz` |

A native component omitted from the file starts at zero. You can edit any
native component in place, for example `fm.Bz .*= 2`.

```julia
fm.Bz                    # longitudinal magnetic field [T]
coord_vec(fm, :z)        # z coordinates [m]
real.(fm.Ez)             # real part of the longitudinal electric field
imag.(fm.Ez)             # imaginary part
abs.(fm.Ez)              # longitudinal component amplitude
field_magnitude(fm, :E)  # electric field magnitude across all three components
field_magnitude(fm, :B)  # magnetic field magnitude across all three components
```

Use [`coord_vecs`](@ref) to retrieve all three coordinate vectors in axis order.

## Grid Information

As in Python `beamphysics`, coordinate vectors start at `gridOriginOffset` and
advance by `gridSpacing`. The reader preserves `gridLowerBound` in `fm.metadata`.
For a rectangular Bmad map, the first sample relative to the element anchor is
`gridOriginOffset + gridLowerBound .* gridSpacing` under Bmad's indexing
convention. When that offset term is nonzero, account for it when comparing
sample locations with Bmad. See the [Bmad manual, §§5.16.4 and 43.2](https://www.classe.cornell.edu/bmad/bmad-manual-2026-03-09.pdf).

```julia
fm.grid_origin     # coordinates of the first grid point
fm.grid_spacing    # grid spacing per axis
fm.grid_size       # number of grid points per axis
maxs(fm)           # grid maximum coordinates
axis_index(fm, :z) # 1-based index of the z axis

# Per-axis helpers (rectangular)
xmin(fm); xmax(fm); dx(fm)
ymin(fm); ymax(fm); dy(fm)
zmin(fm); zmax(fm); dz(fm)

# Per-axis helpers (cylindrical)
rmin(fm); rmax(fm); dr(fm)
zmin(fm); zmax(fm); dz(fm)
```

See [`axis_index`](@ref), [`mins`](@ref), [`maxs`](@ref), [`deltas`](@ref).

## Sampling an axisymmetric solenoid

A cylindrical map with one theta sample represents an axisymmetric field. Use
the stored theta coordinate for interpolation and the particle's azimuth to
rotate the field into Cartesian coordinates. For the LCLS-II solenoid,
the azimuthal magnetic component is zero, so only `Br` and `Bz` enter:

```julia
fm = FieldMesh(joinpath(pkgdir(BeamPhysics), "test", "lcls2_solenoid_fieldmesh.h5"))
br = interpolator(fm, :Br)
bz = interpolator(fm, :Bz)
theta_grid = only(coord_vec(fm, :theta))

x, y, z = 0.012, 0.016, -0.04  # particle position [m]
radius = hypot(x, y)
azimuth = atan(y, x)
radial_field = br(radius, theta_grid, z)
Bx = radial_field * cos(azimuth)
By = radial_field * sin(azimuth)
Bz = bz(radius, theta_grid, z)
```

For a single sample use `interpolate(fm, :Bz, (radius, theta_grid, z))`.
Both [`interpolator`](@ref) and [`interpolate`](@ref) also accept grid arrays.
For a complex RF map, `interpolate(fm, abs.(fm.Ez), point)` interpolates the
amplitude of the grid samples. Use `abs(interpolate(fm, :Ez, point))` for the
amplitude of the interpolated complex field.

## Field Properties

```julia
is_static(fm)        # true if harmonic == 0
is_pure_magnetic(fm) # true if every electric component is zero or absent
is_pure_electric(fm) # true if every magnetic component is zero or absent
frequency(fm)        # oscillation frequency [Hz]
```

See [`is_static`](@ref), [`is_pure_magnetic`](@ref), [`is_pure_electric`](@ref),
[`frequency`](@ref).

## Grid Modification

```julia
set_min!(fm, :z, -0.1)  # shift z origin to -0.1 m
set_max!(fm, :z,  0.1)  # shift z origin so max = 0.1 m
```

See [`set_min!`](@ref), [`set_max!`](@ref).
