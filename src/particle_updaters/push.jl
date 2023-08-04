struct SimpleParticlePush{S,T} <: AbstractSimulationStep
    species::S
    timestep::T
end

function step!(step::SimpleParticlePush)
    species = step.species

    # TODO: this doesn't work for variable weight particles
    # species.positions .= species.positions .+ (step.timestep / species.mass / species.weights[1]) .* species.momentums
    species.positions .=
        species.positions .+ (step.timestep / species.mass) .* species.momentums
    species.momentums .= species.momentums .+ step.timestep .* species.forces
end

struct CommunicateSpecies{S,G} <: AbstractSimulationStep
    species::S
    grid::G
end

function step!(step::CommunicateSpecies)
    for n in eachindex(step.species.positions)
        step.species.positions[n] =
            step.species.positions[n] .-
            floor.(
                (step.species.positions[n] .- step.grid.lower_bounds) ./
                sim_lengths(step.grid)
            ) .* sim_lengths(step.grid)
    end
end
