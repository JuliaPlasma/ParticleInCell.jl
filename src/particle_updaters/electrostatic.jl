"""
    ElectrostaticParticlePush(species, E, timestep, [interpolation_order=1])

An update step that moves and accelerates particles of `species` according
to the electric field `E`. The field will be interpolated to the particle
positions using b-splines of `interpolation_order`.
"""
struct ElectrostaticParticlePush{S,F,T,IF} <: AbstractSimulationStep
    species::S
    E::F
    timestep::T

    interpolation_order::Int
    interpolation_width::Int
    interpolation_function::IF

    function ElectrostaticParticlePush(
        species::S,
        E::F,
        timestep::T,
        interpolation_order = 1,
    ) where {S,F,T}
        interpolation_width = bs_interp_width(interpolation_order)
        interpolation_function = select_bs_interp(interpolation_order)

        new{S,F,T,typeof(interpolation_function)}(
            species,
            E,
            timestep,
            interpolation_order,
            interpolation_width,
            interpolation_function,
        )
    end
end

function step!(step::ElectrostaticParticlePush)
    species = step.species

    for n in eachindex(species)
        # Push the particle based on its current velocity
        particle_position!(
            species,
            n,
            particle_position(species, n) .+
            (step.timestep / particle_mass(species, n)) .* particle_momentum(species, n),
        )

        # Accelerate the particle according to E
        # Find which cell the particle is in, and create a CartesianIndices
        # object that extends +/- interpolation_width in all directions
        particle_cell_coord, Is = phys_coords_to_cell_index_ittr(
            step.E,
            particle_position(species, n),
            step.interpolation_width,
        )

        for I in Is
            grid_cell_coord = cell_index_to_cell_coords(step.E, I, 1)
            dist = Tuple(particle_cell_coord .- grid_cell_coord)
            interp_weights = step.interpolation_function.(dist)
            interp_weight = prod(interp_weights)
            particle_momentum!(
                species,
                n,
                particle_momentum(species, n) .+
                (interp_weight * step.timestep * particle_charge(species, n)) .*
                step.E.values[I],
            )
        end
    end
end
