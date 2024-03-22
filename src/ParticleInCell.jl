module ParticleInCell

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end ParticleInCell

using FFTW
using StaticArrays
using StructEquality
using SparseArrays
using LinearSolve

import Base: eachindex

include("grids/grid.jl")
export AbstractGrid, UniformCartesianGrid

include("species/species.jl")
export AbstractSpecies,
    particle_charge,
    physical_charge,
    particle_mass,
    physical_mass,
    particle_weight,
    particle_position,
    particle_position!,
    particle_momentum,
    particle_momentum!,
    particle_velocity,
    physical_momentum,
    VariableWeightSpecies,
    electrons

include("field.jl")
export Field, NodeOffset, EdgeOffset, FaceOffset, CenterOffset, num_elements

abstract type AbstractSimulationStep end
function step!(::T) where {T<:AbstractSimulationStep}
    error("step! not defined for type $T")
end
export step!

include("field_utils.jl")
export FiniteDifferenceToEdges, AverageEdgesToNodes

include("poisson.jl")
export PoissonSolveFFT

include("linear_solve.jl")
export LinearSolveStep

include("interpolation.jl")
export BSplineChargeInterpolation

include("communicate/communicate_species.jl")
export CommunicateSpecies
include("communicate/communicate_field.jl")
export CommunicateGuardCells

include("particle_updaters/electrostatic.jl")
export ElectrostaticParticlePush
include("particle_updaters/boris.jl")
export BorisParticlePush

include("simulation.jl")
export Simulation, create_electrostatic_simulation

end
