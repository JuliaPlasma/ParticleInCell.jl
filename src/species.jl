abstract type AbstractSpecies end

"""
Stores the positions and momentums of particles that share a common charge and
mass. Each particle can have a different weight.
"""
struct Species{D,T} <: AbstractSpecies
    positions::Vector{SVector{D,T}}
    momentums::Vector{SVector{D,T}}
    forces::Vector{SVector{D,T}}
    weights::Vector{T}
    charge::T
    mass::T
end
