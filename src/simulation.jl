abstract type AbstractSimulationStep end

function step!(sim, step::T) where {T <: AbstractSimulationStep}
    error("step! not defined for type $T")
end

abstract type AbstractSimulationPipeline end

struct SimulationPipeline <: AbstractSimulationPipeline
    steps::Vector{AbstractSimulationStep}
end

struct PeriodicSimulationPipeline <: AbstractSimulationPipeline
    pipeline_period::Int
    steps::Vector{AbstractSimulationStep}
end

abstract type AbstractSimulation end

struct Simulation{T} <: AbstractSimulation
    timestep::T
    time::Val{T}

    fields::Vector{AbstractField}
    species::Vector{AbstractSpecies}

    pipelines::Vector{AbstractSimulationPipeline}
end
