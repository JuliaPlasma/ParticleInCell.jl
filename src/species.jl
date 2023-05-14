abstract type AbstractSpecies end

struct Species{D, T} <: AbstractSpecies
    positions::Vector{SVector{D, T}}
    momentums::Vector{SVector{D, T}}
    forces   ::Vector{SVector{D, T}}
    weights::Vector{T}
    charge::T
    mass::T
end
