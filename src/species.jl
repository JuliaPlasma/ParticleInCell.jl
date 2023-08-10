abstract type AbstractSpecies end

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

"""
    electrons(positions::Matrix, momentums::Matrix, weights::Vector)
    electrons(positions::Matrix, momentums::Matrix, weight)

Creates a `Species` with the given `positions`, `momentums`, and `weight`, and the charge
and mass of a single electron.
"""
function electrons(positions::Matrix{T}, momentums::Matrix{T}, weights::Vector{T}) where {T}
    return Species(positions, momentums, weights, 1.602176634e-19, 9.1093837015e-31)
end

function electrons(positions::Matrix{T}, momentums::Matrix{T}, weight::T) where {T}
    weights = fill(weight, axes(positions, 2))
    return electrons(positions, momentums, weights)
end

function electrons(positions::Vector{T}, momentums::Vector{T}, weights::Vector{T}) where {T}
    position_matrix = reshape(positions, (1, length(positions)))
    momentum_matrix = reshape(momentums, (1, length(momentums)))
    return Species(
        position_matrix,
        momentum_matrix,
        weights,
        1.602176634e-19,
        9.1093837015e-31,
    )
end


function electrons(positions::Vector{T}, momentums::Vector{T}, weight::T) where {T}
    weights = fill(weight, length(positions))
    return electrons(positions, momentums, weights)
end
