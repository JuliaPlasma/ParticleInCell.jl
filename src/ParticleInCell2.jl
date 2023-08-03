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

include("grids/grid.jl")
export AbstractGrid, UniformCartesianGrid, node, edge, face

include("field.jl")
export Field

include("species.jl")
export Field, Species

abstract type AbstractSimulationStep end
function step!(step::T) where {T<:AbstractSimulationStep}
    error("step! not defined for type $T")
end
export step!

include("field_utils.jl")
export FiniteDifferenceToEdges, AverageEdgesToNodes, CommunicateGuardCells

include("poisson.jl")
export PoissonSolveFFT

include("interpolation.jl")
export BSplineChargeInterpolation, BSplineFieldInterpolation

include("push.jl")
export SimpleParticlePush, CommunicateSpecies
include("boris.jl")
export BorisParticlePush

end
