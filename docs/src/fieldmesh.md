# Field Mesh

## Overview

`FieldMesh` stores gridded electromagnetic field data in openPMD-beamphysics format.
It supports both rectangular (`x, y, z`) and cylindrical (`r, theta, z`) geometries,
static and RF fields, and complex field factors.

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

Field components and coordinates are accessed via string indexing:

```julia
fm["Bz"]      # Bz component, scaled by field_scale * exp(i*phase)
fm["Er"]      # Er component (cylindrical only)
fm["z"]       # z coordinate vector
fm["r"]       # r coordinate vector (cylindrical)

# Operator prefixes (useful for complex RF fields)
fm["re_Bz"]   # real part of Bz
fm["im_Bz"]   # imaginary part of Bz
fm["abs_Ez"]  # |Ez|

# Field magnitudes
fm["B"]       # |B| across all components
fm["E"]       # |E| across all components
```

The underlying functions are [`coord_vec`](@ref) and [`scaled_component`](@ref);
use [`coord_vecs`](@ref) to retrieve all three coordinate vectors at once.

See [`coord_vec`](@ref), [`coord_vecs`](@ref), [`scaled_component`](@ref).

## Grid Information

As in Python `beamphysics`, coordinate vectors start at `gridOriginOffset` and
advance by `gridSpacing`. The reader preserves `gridLowerBound` in `fm.metadata`.
For a rectangular Bmad map, the first sample relative to the element anchor is
`gridOriginOffset + gridLowerBound .* gridSpacing` under Bmad's indexing
convention. When that offset term is nonzero, account for it when comparing
sample locations with Bmad. See the [Bmad manual, §§5.16.4 and 43.2](https://www.classe.cornell.edu/bmad/bmad-manual-2026-03-09.pdf).

```julia
mins(fm)            # grid minimum coordinates
maxs(fm)            # grid maximum coordinates
deltas(fm)          # grid spacing per axis
fm.grid_size        # number of grid points per axis
axis_index(fm, "z") # 1-based index of the "z" axis

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
fm = FieldMesh("lcls2_solenoid_fieldmesh.h5")
br = interpolator(fm, "Br")
bz = interpolator(fm, "Bz")
theta_grid = only(fm["theta"])

x, y, z = 0.012, 0.016, -0.04  # particle position [m]
radius = hypot(x, y)
azimuth = atan(y, x)
radial_field = br(radius, theta_grid, z)
Bx = radial_field * cos(azimuth)
By = radial_field * sin(azimuth)
Bz = bz(radius, theta_grid, z)
```

The component keys select fields in the map's native basis. Looking up
`fm["Bx"]` on a cylindrical mesh does not perform this rotation.

## Field Properties

```julia
is_static(fm)        # true if harmonic == 0
is_pure_magnetic(fm) # true if no electric components
is_pure_electric(fm) # true if no magnetic components
scale(fm)            # field_scale factor
phase(fm)            # phase in radians: -2π * RFphase
factor(fm)           # complex factor: scale * exp(i*phase)
frequency(fm)        # oscillation frequency [Hz]
```

See [`is_static`](@ref), [`is_pure_magnetic`](@ref), [`is_pure_electric`](@ref),
[`scale`](@ref), [`phase`](@ref), [`factor`](@ref), [`frequency`](@ref).

## Grid Modification

```julia
set_min!(fm, "z", -0.1)  # shift z origin to -0.1 m
set_max!(fm, "z",  0.1)  # shift z origin so max = 0.1 m
```

See [`set_min!`](@ref), [`set_max!`](@ref).
