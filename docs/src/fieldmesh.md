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

Field properties such as `fm.Ez` give a writable view of the physical samples
in the mesh's native basis. Each view has shape `fm.grid_size` and dimension
order `fm.axis_labels`.

| Mesh geometry | Electric field [V/m] | Magnetic field [T] |
|---|---|---|
| Rectangular | `fm.Ex`, `fm.Ey`, `fm.Ez` | `fm.Bx`, `fm.By`, `fm.Bz` |
| Cylindrical | `fm.Er`, `fm.Etheta`, `fm.Ez` | `fm.Br`, `fm.Btheta`, `fm.Bz` |

A native component omitted from the file starts at zero.

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

## Scaling an RF field map

`fm.field_scale` multiplies every electric and magnetic component.
For a cavity that gives 8 MV at its current setting, multiplying the scale
by 10/8 gives 10 MV at the same phase.

```julia
reference_scale = fm.field_scale
fm.field_scale = reference_scale * (10.0 / 8.0)
```

To change the RF phase, set `fm.rf_phase` in turns. The phasor convention is
`exp(-2π * im * fm.rf_phase)`. For a local profile correction, edit the physical
component directly. For example, `fm.Ez[1, 1, k] *= 0.98` reduces the amplitude
of one sample by 2% at the current scale and phase.

A view of an unphased static map has real elements. If you then apply a phase
to that map, the retained view gives the real part. Access the component again
to get its full complex values.

For a component present in the file, `fm.components[:Ez]` holds the unscaled
profile in its original units, and `fm.component_units[:Ez].scale` gives its
conversion factor to V/m. Changing the declared unit reinterprets the stored
samples immediately, including through previously obtained views.
A physical edit updates the corresponding unscaled
sample. Changing `fm.field_scale` or `fm.rf_phase` leaves the profile intact.
At zero scale all physical samples are zero. To edit the profile then, use
`fm.components[:Ez]` directly or restore a nonzero scale before editing `fm.Ez`.

For example, to reduce the stored longitudinal electric field by 2% across
one transverse plane of an RF cavity, edit that plane in the profile array.
This also works with the cavity switched off.

```julia
k = cld(fm.grid_size[3], 2)     # middle z plane
fm.components[:Ez][:, :, k] .*= 0.98
```

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
Rebuild a reusable interpolator after changing the profile, scale, phase, or grid.
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
