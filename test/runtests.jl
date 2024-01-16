using ParticleInCell
using Test
using ReTestItems

if length(ARGS) == 0
    runtests(ParticleInCell)
else
    runtests(ParticleInCell, tags = Symbol(ARGS[1]))
end
