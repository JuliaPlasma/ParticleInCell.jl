module ParticleInCell2

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    readme_str = read(path, String)

    # Get rid of the title and badges
    index = last(findfirst("\n\n", readme_str)) + 1
    readme_str[index:end]
end ParticleInCell2

using FFTW
using StaticArrays
using StructEquality

include("grids/grid.jl")
export AbstractGrid, UniformCartesianGrid, node, edge, face

include("field.jl")
export Field

include("species.jl")
export Species, electrons

abstract type AbstractSimulationStep end
function step!(step::T) where {T<:AbstractSimulationStep}
    error("step! not defined for type $T")
end
export step!

include("field_utils.jl")
export FiniteDifferenceToEdges, AverageEdgesToNodes

include("poisson.jl")
export PoissonSolveFFT

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

end
