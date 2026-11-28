# [API Reference](@id api)

Complete listing of all exported symbols.

## Types

`Species` is re-exported from AtomicAndPhysicalConstants. For example, `Species("electron")`
selects the electron species when constructing a particle group.

```@docs
ParticleGroup
ParticleGroupView
FieldMesh
PMDUnit
```

## Particle Constructors

```@docs
single_particle
ParticleGroup(h5::Union{HDF5.File,HDF5.Group})
ParticleGroup(file::AbstractString)
```

## Particle Derived Properties

### Group properties

```@docs
nalive
ndead
mass
species_charge
charge
average_current
```

### Per-particle arrays

```@docs
momentum
energy
kinetic_energy
gamma
beta
beta_x
beta_y
beta_z
xp
yp
r
theta
pr
ptheta
Lz
x_bar
px_bar
y_bar
py_bar
Jx
Jy
```

### Emittances

```@docs
norm_emit_x
norm_emit_y
norm_emit_4d
```

## Particle Operations

```@docs
drift!
drift_to_z!
drift_to_t!
set_charge!
in_z_coordinates
in_t_coordinates
```

## Statistics

```@docs
mean(::BeamPhysics.AbstractParticleGroup, ::Any)
std(::BeamPhysics.AbstractParticleGroup, ::Any)
cov(::BeamPhysics.AbstractParticleGroup, ::Any, ::Any)
delta
ptp
norm_emit_calc
twiss_calc
twiss_dispersion
twiss
A_mat_calc
A_inverse_mat_calc
amplitude_calc
particle_amplitude
normalized_particle_coordinate
covariance_matrix
slices
```

## FieldMesh Accessors

```@docs
axis_index
mins
maxs
deltas
coord_vec
coord_vecs
dx
dy
dz
dr
dtheta
xmin
xmax
ymin
ymax
zmin
zmax
rmin
rmax
thetamin
thetamax
frequency
is_static
is_pure_electric
is_pure_magnetic
field_magnitude
interpolator
interpolate
set_min!
set_max!
write_fieldmesh
```

## Units & Display

```@docs
PMDUnit(symbol::String)
nice_scale_prefix
nice_array
limits
plottable_array
pg_units
mathlabel
texlabel
```

## I/O

```@docs
particle_paths
field_paths
load_field_attrs
is_constant_component
constant_component_value
component_data
offset_component_name
particle_array
all_components
component_str
pmd_init
pmd_field_init
write_component_data
write_pmd_bunch
write_pmd_field
particle_data_dict
write_particle_group
```

## Plotting

```@docs
density_plot
marginal_plot
slice_plot
density_and_slice_plot
fieldmesh_plot
fieldmesh_plot_onaxis
```
