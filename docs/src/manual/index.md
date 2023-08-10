# Manual

## Grids
```@docs
AbstractGrid
UniformCartesianGrid
```
### Grid utility functions
```@docs
ParticleInCell2.cell_lengths
ParticleInCell2.sim_lengths
ParticleInCell2.cell_coords_to_phys_coords
ParticleInCell2.phys_coords_to_cell_coords
```

## Species
```@docs
AbstractSpecies
VariableWeightSpecies
electrons
```

### Species utility functions
```@docs
particle_charge
physical_charge
particle_mass
physical_mass
particle_weight
particle_position
particle_position!
particle_momentum
particle_momentum!
particle_velocity
physical_momentum
eachindex
```

## Fields
```@docs
Field
```

## Update steps
```
ParticleInCell2.AbstractUpdateStep
step!
```
