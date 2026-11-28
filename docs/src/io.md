# I/O

## File Format

This package reads and writes files in the
[openPMD-beamphysics](https://github.com/openPMD/openPMD-standard/blob/upcoming-2.0.0/EXT_BeamPhysics.md) HDF5 format.

## Reading Particles

```julia
using BeamPhysics

# Load particles (file must contain exactly one particle path)
pg = ParticleGroup("particles.h5")

# Inspect available paths in a file
using HDF5
h5open("particles.h5", "r") do h5
    println(particle_paths(h5))
end
```

Each entry in `pg.weight` is the macroparticle charge in coulombs.
Particle IDs are available in `pg.id` and optional columns in `pg.extra`,
for example `pg.extra[:p0c]` in eV/c and `pg.extra[:spin_x]` in units of hbar.
The reader adds time offsets by default. Pass `include_time_offset=false` to read the
stored time without its offset.

## Writing Particles

```julia
write_particle_group("output.h5", pg)
# If pg.t holds raw times, pass t_offset in seconds to write a separate timeOffset record.
```

[`write_particle_group`](@ref) preserves particle IDs and the following optional
columns in `pg.extra` when present: `s`, `ix_branch`, `ix_ele`, `location`,
`p0c`, `spin_x`, `spin_y`, `spin_z`, and `charge_state`. Other entries,
such as a user-defined `pg.extra[:diagnostic]`, are in-memory only.

## Reading Field Meshes

```julia
fm = FieldMesh("fieldmesh.h5")

# Inspect available paths
h5open("fieldmesh.h5", "r") do h5
    println(field_paths(h5))
end
```

See [`FieldMesh`](@ref), [`field_paths`](@ref).

Fields load in (x, y, z) or (r, theta, z) order, including files with a reversed
axis layout or the legacy Fortran order marker. Static fields are real; RF fields
retain their complex values. Additional field attributes are retained in `fm.metadata`.

## Writing Field Meshes

```julia
write_fieldmesh("output.h5", fm)
```

See [`write_fieldmesh`](@ref).

## Low-Level HDF5 API

For advanced use, the following functions are also exported. See the
[API Reference](@ref api) for full signatures.

**Reading:** [`particle_paths`](@ref), [`field_paths`](@ref),
[`load_field_attrs`](@ref), [`is_constant_component`](@ref),
[`constant_component_value`](@ref), [`component_data`](@ref),
[`offset_component_name`](@ref), [`particle_array`](@ref),
[`all_components`](@ref), [`component_str`](@ref)

**Writing:** [`pmd_init`](@ref), [`pmd_field_init`](@ref),
[`write_component_data`](@ref), [`write_pmd_bunch`](@ref),
[`write_pmd_field`](@ref), [`particle_data_dict`](@ref),
[`write_particle_group`](@ref)
