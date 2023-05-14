abstract type AbstractSpecies end

# TODO: switch to Vector{NTuple{D, T}}?

struct Species{T} <: AbstractSpecies
    positions::Array{T, 2}
    momentums::Array{T, 2}
    forces::Array{T, 2}
    weights::Array{T, 2}
    charge::T
    mass::T
end
