# BeamPhysics.jl

BeamPhysics.jl is a Julia package for analyzing charged-particle beams and electromagnetic field maps in the openPMD standard, beamphysics extension.

[openPMD BeamPhysics extension specification](https://github.com/openPMD/openPMD-standard/blob/upcoming-2.0.0/EXT_BeamPhysics.md)

## Features

- Load and save particle data (`ParticleGroup`) from openPMD HDF5 files
- Load and save gridded electromagnetic field data (`FieldMesh`) from openPMD HDF5 files
- Compute derived particle properties (emittance, Twiss parameters, etc.)
- Drift particles in free space
- Optional plotting via [Plots.jl](https://github.com/JuliaPlots/Plots.jl)

## Installation

```julia
using Pkg
Pkg.add("BeamPhysics")
```

## Quick Start

```julia
using BeamPhysics

# Load particles from an openPMD HDF5 file
pg = ParticleGroup(joinpath(pkgdir(BeamPhysics), "test", "initial_particles.h5"))

# Access particle properties
pg.x           # x positions [m]
std(pg, :x)     # std of x [m]
norm_emit_x(pg) # normalized emittance [m]

# Load a field mesh
fm = FieldMesh(joinpath(pkgdir(BeamPhysics), "test", "lcls2_solenoid_fieldmesh.h5"))
fm.Bz             # longitudinal magnetic field on the grid [T]
coord_vec(fm, :z) # z coordinate vector [m]
```

## Contents

```@contents
Pages = ["particles.md", "fieldmesh.md", "statistics.md", "io.md", "plotting.md", "api.md"]
Depth = 2
```
