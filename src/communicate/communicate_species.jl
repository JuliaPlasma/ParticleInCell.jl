struct CommunicateSpecies{S,G} <: AbstractSimulationStep
    species::S
    grid::G
end

# TODO...
function step!(step::CommunicateSpecies)
    species = step.species
    for n in eachindex(species)
        particle_position!(
            species,
            n,
            particle_position(species, n) .-
            floor.(
                (particle_position(species, n) .- step.grid.lower_bounds) ./
                sim_lengths(step.grid)
            ) .* sim_lengths(step.grid),
        )
    end
end
