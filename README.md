# ParticleInCell2.jl
<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://USER_NAME.github.io/PACKAGE_NAME.jl/stable) -->
[![Dev docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://adamslc.github.io/ParticleInCell2.jl/dev)
[![CI](https://github.com/adamslc/ParticleInCell2.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/adamslc/ParticleInCell2.jl/actions/workflows/CI.yml)
[![codecov.io](http://codecov.io/github/adamslc/ParticleInCell2.jl/coverage.svg?branch=main)](http://codecov.io/github/adamslc/ParticleInCell2.jl?branch=main)

`ParticleInCell2.jl` is a Julia package for kinetic plasma physics simulation.
Currently, this package is in in a pre-1.0.0 state, and thus breaking changes
should be expected. However, this also means that I am willing to entertain
radical suggestions to improve the functionality of the package. If you are
interested in using `ParticleInCell2` for your plasma research, and you find
that it does not meet you needs, please reach out on either GitHub, or over
email, so that we can discuss how the package can be modified to suite your
needs.

## Getting Started
`ParticleInCell2` is currently not registered in the Julia package registry.
Thus, to install this package, you should `dev` it:
```julia
using Pkg
Pkg.develop(url="https://github.com/adamslc/ParticleInCell2.jl")
```

## Documentation
You can view the latest documentation
[here](https://adamslc.github.io/ParticleInCell2.jl/dev).

## Goals
 * Fast: aim to have core time of less than 1 microsecond per particle per step
   without collisions.
 * Flexible: it should be possible to implement essentially any kinetic plasma
   simulation in `ParticleInCell2.jl`. For common types of simulations, this
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

## Todo
 * New examples:
   * Two-stream instability (check frequency)
   * Capacitively coupled plasma
   * EM radiation from a dipole antenna
 * Add `LinearFieldSolve` step
   * Figure out a flexible way to define a stencil
 * Get rid of force field in Species, and incorporate the field interpolation into
   the particle push?
   * Add Boris push
 * Add electromagnetic update step
 * Add Monte-Carlo collisions
