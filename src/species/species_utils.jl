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
