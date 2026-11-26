# Particles

## Overview

`ParticleGroup` is the primary data structure for beam particles. It stores particle
coordinates, momenta, time, status, weight, and species. Array fields are parameterized
so they can be `Vector`, `CuArray`, or any `AbstractVector`, enabling GPU computation.

See [`ParticleGroup`](@ref) in the [API Reference](@ref api) for the full type definition.

## Creating a ParticleGroup

```julia
using BeamPhysics

# Load from an openPMD HDF5 file
pg = read_particle_group("particles.h5")

# Create a single test particle
pg = single_particle(z=0.01, pz=1e6)

# Construct directly
using AtomicAndPhysicalConstants: Species
pg = ParticleGroup(; x, px, y, py, z, pz, t, status, weight, species=Species("electron"))
```

See [`single_particle`](@ref) and [`read_particle_group`](@ref).

## Accessing Properties

Stored columns are available as properties (for example, `pg.x` and `pg.weight`).
Derived quantities are functions: `energy(pg)` returns one value per particle,
while `charge(pg)` and `norm_emit_x(pg)` return one value for the group.

See the [statistics guide](statistics.md) for weighted group reductions.

## Selecting Particles

```julia
pg[1:100]        # first 100 particles (returns a copy)
view(pg, 1:100)  # non-owning view into the first 100 particles
```

See [`ParticleGroupView`](@ref).

## Joining

```julia
pg3 = pg1 + pg2  # concatenate two ParticleGroups of the same species
```

## In-place Operations

```julia
drift!(pg, 1e-9)       # drift by 1 ns
drift_to_z!(pg, 0.5)   # drift all particles to z = 0.5 m
drift_to_t!(pg)        # drift all particles to their mean time
set_charge!(pg, 1e-12) # rescale weights to total charge = 1 pC
```

## Splitting

```julia
groups = slices(pg; by=:z, n=100)
```

## Derived Properties

Use the following functions to calculate derived quantities.

| Function call | Returns | Units |
|---|---|---|
| `nalive(pg)` | number of alive particles | — |
| `ndead(pg)` | number of dead particles | — |
| `charge(pg)` | total charge | C |
| `mass(pg)` | particle rest mass | eV |
| `species_charge(pg)` | elementary charge of species | C |
| `average_current(pg)` | average current | A |
| `momentum(pg)` | total momentum per particle | eV/c |
| `energy(pg)` | total energy per particle | eV |
| `kinetic_energy(pg)` | kinetic energy per particle | eV |
| `gamma(pg)` | Lorentz factor per particle | — |
| `beta(pg)` | `\|v\|/c` per particle | — |
| `beta_x(pg)` | `vx/c` per particle | — |
| `beta_y(pg)` | `vy/c` per particle | — |
| `beta_z(pg)` | `vz/c` per particle | — |
| `xp(pg)` | `px/pz` per particle | rad |
| `yp(pg)` | `py/pz` per particle | rad |
| `r(pg)` | cylindrical radius per particle | m |
| `theta(pg)` | azimuthal angle per particle | rad |
| `pr(pg)` | radial momentum per particle | eV/c |
| `ptheta(pg)` | azimuthal momentum per particle | eV/c |
| `Lz(pg)` | angular momentum per particle | m·eV/c |
| `norm_emit_x(pg)` | normalized x emittance | m |
| `norm_emit_y(pg)` | normalized y emittance | m |
| `norm_emit_4d(pg)` | normalized 4D emittance | m² |
| `x_bar(pg)` | normalized x coordinate per particle | m^(1/2) |
| `px_bar(pg)` | normalized px coordinate per particle | m^(1/2) |
| `y_bar(pg)` | normalized y coordinate per particle | m^(1/2) |
| `py_bar(pg)` | normalized py coordinate per particle | m^(1/2) |
| `Jx(pg)` | x action per particle | m |
| `Jy(pg)` | y action per particle | m |
