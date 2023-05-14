module ParticleInCell2
    using FFTW
    using StaticArrays

    include("grid.jl")
    include("field.jl")
    include("species.jl")
    export UniformCartesianGrid, Field, Species

    include("simulation.jl")
    export step!

    include("poisson.jl")
    export PoissonSolveFFT

    include("scatter.jl")
    export BSplineChargeInterpolation

    include("gather.jl")
    export BSplineFieldInterpolation

    include("field_utils.jl")
    export FiniteDifferenceToEdges, AverageEdgesToNodes, CommunicateGuardCells
end
