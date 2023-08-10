"""
    VariableWeightSpecies(positions::Matrix, momentums::Matrix, weight::Vector charge, mass)

Stores the positions and momentums of particles that share a common charge and mass. Each
particle can have a different weight, and the momentum space can have a larger dimension
than the position space.
"""
@struct_hash_equal struct VariableWeightSpecies{D,V,T} <: AbstractSpecies{D,V,T}
    positions::Vector{SVector{D,T}}
    momentums::Vector{SVector{V,T}}
    weights::Vector{T}
    charge::T
    mass::T
end

function VariableWeightSpecies(
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

    return VariableWeightSpecies{D,V,T}(pos_vec, mom_vec, weights, charge, mass)
end

particle_charge(species::VariableWeightSpecies, idx) = species.charge * species.weights[idx]
physical_charge(species::VariableWeightSpecies, idx) = species.charge
particle_mass(species::VariableWeightSpecies, idx) = species.mass * species.weights[idx]
physical_mass(species::VariableWeightSpecies, idx) = species.mass
particle_weight(species::VariableWeightSpecies, idx) = species.weights[idx]

particle_position(species::VariableWeightSpecies, idx) = species.positions[idx]
particle_position!(species::VariableWeightSpecies, idx, value) =
    species.positions[idx] = value
particle_momentum(species::VariableWeightSpecies, idx) = species.momentums[idx]
particle_momentum!(species::VariableWeightSpecies, idx, value) =
    species.momentums[idx] = value
# TODO: this is not relativistic...
particle_velocity(species::VariableWeightSpecies, idx) =
    species.momentums[idx] / particle_mass(species, idx)
physical_momentum(species::VariableWeightSpecies, idx) =
    species.momentums[idx] / species.weights[idx]
eachindex(species::VariableWeightSpecies) = eachindex(species.positions)
