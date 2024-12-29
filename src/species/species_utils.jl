"""
    electrons(positions::Matrix, momentums::Matrix, weights::Vector)
    electrons(positions::Matrix, momentums::Matrix, weight)

Creates an `AbsractSpecies` with the given `positions`, `momentums`, and `weight`, and the charge
and mass of a single electron.
"""
function electrons(positions::Matrix{T}, momentums::Matrix{T}, weights::Vector{T}) where {T}
    return VariableWeightSpecies(
        positions,
        momentums,
        weights,
        1.602176634e-19,
        9.1093837015e-31,
    )
end

function electrons(positions::Matrix{T}, momentums::Matrix{T}, weight::T) where {T}
    weights = fill(weight, axes(positions, 2))
    return electrons(positions, momentums, weights)
end

function electrons(positions::Vector{T}, momentums::Vector{T}, weights::Vector{T}) where {T}
    position_matrix = reshape(positions, (1, length(positions)))
    momentum_matrix = reshape(momentums, (1, length(momentums)))
    return VariableWeightSpecies(
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

function quiet_velocities(num_macroparticles, distribution, T=Float64)
    velocities = Vector{T}(undef, num_macroparticles)
    for i in 1:num_macroparticles
        velocities[i] = quantile(distribution, (i - 0.5) / num_macroparticles)
    end
    return velocities
end

function bit_reversed_sequence(n, base=2, T=Float64)
    seq = Vector{T}(undef, n)

    m = round(Int, log(base, n), RoundUp)

    for i in 0:n-1
        j = 0
        for (k, digit) in enumerate(digits(i, base=base, pad=m))
            j += digit * base^(m - k)
        end

        seq[i+1] = j / base^m
    end

    return seq
end
