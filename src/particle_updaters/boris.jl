"""
    BorisParticlePush(species, E, B, timestep)

This update step moves the particles of `species` subject to both an electric
field, `E`, and a magnetic field, `B`. The method is frequently referred to by
the shorthand accelerate-rotate-accelerate because the acceleration from the
electric field is split in half, and applied before and after the magnetic
field rotation.

An applied magnetic field forces charged particles to travel in circular orbits
perpendicular to the magnetic field. Thus, a simulation using the Boris method
only makes sense when there are at least two velocity components, and in this
case, the applied magnetic field must be strictly perpendicular to the velocity
components. So for example, one could have a 1d2v simulation with a spatial
grid along the x-axis, and velocity components in the x and y directions. In
this case, the magnetic field would be forced to point along the z axis.
If all three velocity components are included in the simulation, then the
magnetic field can point in any direction.

For more details on the method, see [Sections 4.3 and 4.4 of Birdsall and
Langdon](@cite birdsall2004), or the [Boris' original conference
proceedings](@cite boris1970).
"""
struct BorisParticlePush{V,S,F,T,IF} <: AbstractSimulationStep
    species::S
    E::F
    B::F
    timestep::T

    interpolation_order::Int
    interpolation_width::Int
    interpolation_function::IF

    function BorisParticlePush(
        species::Species{D,V,T},
        E::F,
        B::F,
        timestep::TS,
        interpolation_order = 1,
    ) where {D,V,T,F,TS}
        interpolation_width = bs_interp_width(interpolation_order)
        interpolation_function = select_bs_interp(interpolation_order)

        new{V,Species{D,V,T},F,TS,typeof(interpolation_function)}(
            species,
            E,
            B,
            timestep,
            interpolation_order,
            interpolation_width,
            interpolation_function,
        )
    end
end

function step!(step::BorisParticlePush{1})
    error("Cannot apply Boris particle update to a species with one velocity dimension.")
end

function step!(step::BorisParticlePush{2})
    species = step.species
    # TODO: this will break for an empty species...
    pos_dim = length(species.positions[1])

    for n in eachindex(species.positions)
        # Push the particle based on its current velocity
        species.positions[n] =
            species.positions[n] .+
            (step.timestep / species.mass / species.weights[n]) .*
            species.momentums[n][1:pos_dim]

        # Accelerate the particle according to E and B
        # The Boris algorithm first does half the acceleration by E, then the
        # entire v x B rotation, and then the second half of the acceleration.

        # Find which cell the particle is in, and create a CartesianIndices
        # object that extends +/- interpolation_width in all directions
        particle_cell_coord, Is = phys_coords_to_cell_index_ittr(
            step.E,
            species.positions[n],
            step.interpolation_width,
        )

        # Calculate the interpolated electric and magnetic fields at the
        # location of the particle
        # TODO: this won't work if E and B are located at different places on
        # the grid cell...
        local_E_partial = zeros(axes(step.E.values)[end])
        local_B = 0.0 # B must be a 1D field
        for I in Is
            grid_cell_coord = cell_index_to_cell_coords(step.E, I)
            dist = Tuple(particle_cell_coord .- grid_cell_coord)
            interp_weights = step.interpolation_function.(dist)
            interp_weight = prod(interp_weights)
            local_E_partial .+= interp_weight .* step.E.values[I]
            local_B += interp_weight * step.B.values[I]
        end

        # TODO: this can probably be optimized
        # The spatial dimension must be either 1 or 2. If it is 1, we need
        # to pad the electric field with a zero, because there is no field
        # in the extra velocity direction.
        if length(local_E_partial) < 2
            local_E = (local_E_partial..., 0.0)
        else
            local_E = local_E_partial
        end

        # Apply the first half of the acceleration
        species.momentums[n] =
            species.momentums[n] .+
            (step.timestep * species.charge * species.weights[n] / 2) .* local_E

        # With two velocity dimensions, the magnetic field must be
        # perpendicular to both fields, and thus we can use a shortcut due to
        # Buneman to calculate the rotation. See section 4.3 of Birdsall and
        # Langdon for details. I use the same notation here.
        t = species.charge * local_B * step.timestep / (species.mass * 2)
        s = 2 * t / (1 + t^2)

        px_minus = species.momentums[n][1]
        py_minus = species.momentums[n][2]
        px_prime = px_minus + py_minus * t
        py_plus = py_minus - px_prime * s
        px_plus = px_prime + py_plus * t

        species.momentums[n] = SVector(px_plus, py_plus)

        # Apply the second half of the acceleration
        species.momentums[n] =
            species.momentums[n] .+
            (step.timestep * species.charge * species.weights[n] / 2) .* local_E
    end
end

function step!(step::BorisParticlePush{3})
    error("Not yet implemented...")
end
