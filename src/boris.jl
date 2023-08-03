"""
    BorisParticlePush(species, E, B, timestep)

This update step moves the particles of `species` subject to both an electric
field, `E`, and a magnetic field, `B`. The method is frequently referred to by
the shorthand accelerate-rotate-accelerate because the acceleration from the
electric field is split in half, and applied before and after the magnetic
field rotation.

For more details on the method, see [birdsall2004](@citet).
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
