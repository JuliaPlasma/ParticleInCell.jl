abstract type AbstractSpecies end

"""
Stores the positions and momentums of particles that share a common charge and
mass. Each particle can have a different weight.
"""
struct Species{D,V,T} <: AbstractSpecies
    positions::Vector{SVector{D,T}}
    momentums::Vector{SVector{V,T}}
    weights::Vector{T}
    charge::T
    mass::T
end
