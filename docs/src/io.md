# I/O

## File Format

This package reads and writes files in the
[openPMD-beamphysics](https://github.com/openPMD/openPMD-standard/blob/upcoming-2.0.0/EXT_BeamPhysics.md) HDF5 format.

## Reading Particles

```julia
using BeamPhysics

# Load a file with one particle path and one species
pg = ParticleGroup("particles.h5")
```

For a file with multiple iterations, list the particle paths and select the
one you want. For example, if the list includes `"/data/1/"`:

```julia
using HDF5

paths = h5open(particle_paths, "particles.h5", "r")
pg = h5open("particles.h5", "r") do h5
    ParticleGroup(h5["/data/1/"])
end
```

If the selected group contains multiple species, inspect its children with
`keys(h5["/data/1/"])` inside the block and select the species explicitly,
for example `ParticleGroup(h5["/data/1/electron"])`.

Each entry in `pg.weight` is the macroparticle charge in coulombs.
Particle IDs are available in `pg.id` and optional columns in `pg.extra`,
for example `pg.extra[:p0c]` in eV/c and `pg.extra[:spin_x]` in units of hbar.
The reader adds time offsets by default. Pass `include_time_offset=false` to read the
stored time without its offset.

## Writing Particles

```julia
write_particle_group("particles_out.h5", pg)
# If pg.t holds raw times, pass t_offset in seconds to write a separate timeOffset record.
```

[`write_particle_group`](@ref) preserves particle IDs and the following optional
columns in `pg.extra` when present: `s`, `ix_branch`, `ix_ele`, `location`,
`p0c`, `spin_x`, `spin_y`, `spin_z`, and `charge_state`. Other entries,
such as a user-defined `pg.extra[:diagnostic]`, are in-memory only.

## Reading Field Meshes

```julia
using BeamPhysics

fm = FieldMesh("fieldmesh.h5")  # file with one field map
```

For a file with multiple field maps, list their paths and load the group you
want. For example, if the list includes `"/ExternalFieldMesh/1/"`:

```julia
using HDF5

paths = h5open(field_paths, "fieldmesh.h5", "r")
fm = h5open("fieldmesh.h5", "r") do h5
    FieldMesh(h5["/ExternalFieldMesh/1/"])
end
```

See [`FieldMesh`](@ref), [`field_paths`](@ref).

Fields load in (x, y, z) or (r, theta, z) order, including files with a reversed
axis layout or the legacy Fortran order marker. The reader retains the unscaled
samples in `fm.components` and each component's openPMD unit in `fm.component_units`.
Each unit has its own `unitSI` conversion factor. The mesh-wide `fieldScale`
and `RFphase` attributes are available as `fm.field_scale` and `fm.rf_phase`.
Physical properties such as `fm.Bz` apply the component unit, scale, and phase
on access. Additional file attributes are kept in `fm.metadata`.

## Writing Field Meshes

```julia
write_fieldmesh("fields_out.h5", fm)
```

[`write_fieldmesh`](@ref) writes each component's unscaled samples and its own unit
attributes, together with the current mesh scale and phase. A map written with
zero `fieldScale` retains its underlying profile when reloaded.

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
