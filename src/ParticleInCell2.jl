module ParticleInCell2
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
end
