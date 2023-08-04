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
