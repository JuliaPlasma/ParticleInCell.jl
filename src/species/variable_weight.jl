"""
    Species(positions::Matrix, momentums::Matrix, weight::Vector charge, mass)

Stores the positions and momentums of particles that share a common charge and mass. Each
particle can have a different weight, and the momentum space can have a larger dimension
than the position space.
"""
@struct_hash_equal struct Species{D,V,T} <: AbstractSpecies
    positions::Vector{SVector{D,T}}
    momentums::Vector{SVector{V,T}}
    weights::Vector{T}
    charge::T
    mass::T
end

function Species(
    positions::Matrix{T},
    momentums::Matrix{T},
    weights::Vector{T},
    charge::T,
    mass::T,
) where {T}
    D, pos_nparticles = length.(axes(positions))
    V, mom_nparticles = length.(axes(momentums))

    @assert pos_nparticles == mom_nparticles == length(weights)

    pos_vec = Vector{SVector{D,T}}(undef, pos_nparticles)
    mom_vec = Vector{SVector{V,T}}(undef, pos_nparticles)
    for i = 1:pos_nparticles
        pos_vec[i] = SVector(positions[:, i]...)
        mom_vec[i] = SVector(momentums[:, i]...)
    end

    return Species{D,V,T}(pos_vec, mom_vec, weights, charge, mass)
end

particle_charge(species::Species, idx) = species.charge * species.weights[idx]
physical_charge(species::Species, idx) = species.charge
particle_mass(species::Species, idx) = species.mass * species.weights[idx]
physical_mass(species::Species, idx) = species.mass
particle_weight(species::Species, idx) = species.weights[idx]

particle_position(species::Species, idx) = species.positions[idx]
particle_position!(species::Species, idx, value) = species.positions[idx] = value
particle_momentum(species::Species, idx) = species.momentums[idx]
particle_momentum!(species::Species, idx, value) = species.momentums[idx] = value
# TODO: this is not relativistic...
particle_velocity(species::Species, idx) =
    species.momentums[idx] / particle_mass(species, idx)
physical_momentum(species::Species, idx) = species.momentums[idx] / particle.weights[idx]
eachindex(species::Species) = eachindex(species.positions)