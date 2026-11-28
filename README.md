# BeamPhysics.jl

[![Read the documentation](https://img.shields.io/badge/Read_the_documentation-4063D8?style=for-the-badge)](https://ndwang.github.io/BeamPhysics.jl/)

BeamPhysics.jl is a Julia package for analyzing charged-particle beams and electromagnetic field maps in the openPMD standard, beamphysics extension.

[openPMD BeamPhysics extension specification](https://github.com/openPMD/openPMD-standard/blob/upcoming-2.0.0/EXT_BeamPhysics.md)

## Installation

```julia
using Pkg
Pkg.add("BeamPhysics")
```

## Quick Start

```julia
using BeamPhysics

pg = ParticleGroup("path/to/particles.h5")
std(pg, :x)      # charge-weighted RMS beam size in x [m]
norm_emit_x(pg)  # normalized horizontal emittance [m]

fm = FieldMesh("path/to/fieldmap.h5")
fm.Bz             # longitudinal magnetic field on the grid [T]
```
