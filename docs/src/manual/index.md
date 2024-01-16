# Manual

## Grids
```@docs
AbstractGrid
UniformCartesianGrid
```
### Grid utility functions
```@docs
ParticleInCell.cell_lengths
ParticleInCell.sim_lengths
ParticleInCell.cell_coords_to_phys_coords
ParticleInCell.phys_coords_to_cell_coords
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
ParticleInCell.AbstractUpdateStep
step!
```

```@docs
ElectrostaticParticlePush
BorisParticlePush
```

## Misc
```@docs
ParticleInCell.compute_knots
ParticleInCell.compute_bspline_coeffs
```
