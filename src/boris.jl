"""
    BorisParticlePush(species, E, B, timestep)

This update step moves the particles of `species` subject to both an electric
field, `E`, and a magnetic field, `B`. The method is frequently referred to by
the shorthand accelerate-rotate-accelerate because the acceleration from the
electric field is split in half, and applied before and after the magnetic
field rotation.

For more details on the method, see [Sections 4.3 and 4.4 of Birdsall and
Langdon](@cite birdsall2004), or the [Boris' original conference
proceedings](@cite boris1970).
"""
struct BorisParticlePush{S, F, T} <: AbstractSimulationStep
    species::S
    E::F
    B::F
    timestep::T
end

function step!(step::BorisParticlePush)
    return
end
