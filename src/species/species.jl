"""
Subtypes of `AbstractSpecies` represent a group of macroparticles with a single physical
charge and mass. Each particle of a species can optionally represent a different number of
physical particles.
"""
abstract type AbstractSpecies{D,V,T} end

"""
    particle_charge(species, idx)

Returns the charge of the macroparticle of `species` with index `idx`. See also
[`physical_charge`](@ref).
"""
function particle_charge end

"""
    physical_charge(species, idx)

Returns the charge of a physical particle represented by the macroparticle of `species` with
index `idx`. See also [`particle_charge`](@ref).
"""
function physical_charge end

"""
    particle_mass(species, idx)

Returns the mass of the macroparticle of `species` with index `idx`. See also
[`physical_mass`](@ref).
"""
function particle_mass end

"""
    physical_mass(species, idx)

Returns the mass of a physical particle represented by the macroparticle of `species` with
index `idx`. See also [`particle_charge`](@ref).
"""
function physical_mass end

"""
    particle_weight(species, idx)

Returns the number of  physical particles represented by the macroparticle of `species` with
index `idx`.
"""
function particle_weight end

"""
    particle_position(species, idx)

Returns the position of the macroparticle of `species` with index `idx`. See also
[`particle_momentum`](@ref) and [`particle_position!`](@ref).
"""
function particle_position end

"""
    particle_position!(species, idx, value)

Sets the position of the macroparticle of `species` with index `idx` to `value`. See also
[`particle_momentum!`](@ref) and [`particle_position`](@ref).
"""
function particle_position! end

"""
    particle_momentum(species, idx)

Returns the momentum of the macroparticle of `species` with index `idx`. See also
[`particle_momentum!`](@ref), [`particle_velocity`](@ref), and [`physical_momentum`](@ref).
"""
function particle_momentum end

"""
    particle_momentum!(species, idx, value)

Sets the momentum of the macroparticle of `species` with index `idx` to `value`. See also
[`particle_momentum`](@ref), [`particle_velocity`](@ref), and [`physical_momentum`](@ref).
"""
function particle_momentum! end

"""
    particle_velocity(species, idx)

Returns the velocity of the macroparticle of `species` with index `idx`. See also
[`particle_momentum`](@ref), [`particle_momentum!`](@ref), and [`physical_momentum`](@ref).
"""
function particle_velocity end

"""
    physical_momentum(species, idx)

Returns the momentum of a physical particle represented by the macroparticle of `species`
with index `idx`. See also [`particle_momentum`](@ref), [`particle_momentum!`](@ref), and
[`particle_velocity`](@ref).
"""
function physical_momentum end

"""
    eachindex(species)

Returns an iterator to the index of each of the macroparticles contained in `species`.
"""
function eachindex end

include("variable_weight.jl")
include("species_utils.jl")
