module ParticleInCell2
using FFTW
using StaticArrays

include("grid.jl")
include("field.jl")
include("species.jl")
export UniformCartesianGrid, Field, Species

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
end
