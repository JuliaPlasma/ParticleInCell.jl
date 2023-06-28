# ParticleInCell2.jl: Kinetic Plasma Physics in Julia

`ParticleInCell2` is a package for doing plasma physics in Julia. Specifically,
it focuses on the simulation of kinetic (non-thermal) plasmas using
particle-in-cell (PIC) algorithms.

## Installation
`ParticleInCell2` is currently not registered in the Julia package registry.
Thus, to install this package, you should `dev` it:
```julia
using Pkg
Pkg.develop(url="https://github.com/adamslc/ParticleInCell2.jl")
```

Once installed, you can load the package using
```julia
using ParticleInCell2
```

## Documentation sections
The `ParticleInCell2` documentation is divided into the following sections:
- [Tutorials](@ref): step-by-step guides for basic package functionality.
- [How-to Guides](@ref): more complex examples of package capabilities.
- [Discussions](@ref): overviews of plasma simulation in general, and the specific
  design of the `ParticleInCell2` package.
- [Reference Manual](@ref): detailed description of each public function and
  type in the package.
