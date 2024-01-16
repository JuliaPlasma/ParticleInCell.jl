![ParticleInCell.jl Logo](https://raw.githubusercontent.com/JuliaPlasma/ParticleInCell.jl/main/logo/logo.gif)

[![Latest Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaplasma.github.io/ParticleInCell.jl/dev)
[![CI Status](https://github.com/JuliaPlasma/ParticleInCell.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaPlasma/ParticleInCell.jl/actions/workflows/CI.yml)
[![Code Coverage Statistics](https://codecov.io/gh/adamslc/ParticleInCell2.jl/branch/main/graph/badge.svg)](http://codecov.io/github/adamslc/ParticleInCell2.jl)

`ParticleInCell.jl` is a Julia package for kinetic plasma physics simulation.
Specifically, it focuses on the simulation of kinetic (non-thermal) plasmas
using particle-in-cell (PIC) algorithms. Currently, this package is in in a
pre-1.0.0 state, and thus breaking changes should be expected. However, this
also means that I am willing to entertain radical suggestions to improve the
functionality of the package. If you are interested in using `ParticleInCell`
for your plasma research, and you find that it does not meet you needs, please
reach out on either GitHub, or over email, so that we can discuss how the
package can be modified to suite your needs.

## Getting Started
`ParticleInCell` is currently not registered in the Julia package registry.
Thus, to install this package, you should use `Pkg.develop`:
```julia
using Pkg
Pkg.develop(url="https://github.com/JuliaPlasma/ParticleInCell.jl")
```

## Documentation
You can view the latest documentation
[here](https://adamslc.github.io/ParticleInCell2.jl/dev).

## Goals
 * Fast: aim to have core time of less than 1 microsecond per particle per step
   without collisions.
 * Flexible: it should be possible to implement essentially any kinetic plasma
   simulation in `ParticleInCell.jl`. For common types of simulations, this
   might mean just piecing together components that are already included. More
   esoteric problems might require writing custom types that implement the
   desired algorithms. The advantage of writing this package in Julia is that
   these custom types will be just as performant as native components that are
   included in the package.
 * Scalable: the eventual goal is to enable scaling across an essentially
   unlimited number of cores using Julia's native multithreading for
   parallelization on a single node, and `MPI.jl` for communication across
   nodes. The goal is to support two different modes of parallelization:
   * Each core is responsible for a single rectangular subdomain. The domain
     assigned to an entire node is also rectangular, which imposes constraints
     on how the node domain can be subdivided into subdomains for each core.
     Load balancing is achieved by varying the relative sizes of the domains
     such that each core has a similar amount of work per step.
   * The simulation domain is subdivided into subdomains called 'patches', and
     every node is assigned a list of patches that it is responsible for
     updating. The cores on each node work collaboratively on the list, each
     choosing one patch to work on, and then selecting another when they are
     finished. Load balancing is achieved by swapping patches between nodes to
     balance the workload while also seeking to minimize communication time by
     keeping the surface area of each node's responsibilities minimized. In
     order for this scheme to effectively load balance, it must be the case that
     the total number of patches is larger (ideally much larger) that the total
     number of cores.
