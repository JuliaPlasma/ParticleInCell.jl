using ParticleInCell2
using Test
using ReTestItems

if length(ARGS) == 0
    runtests(ParticleInCell2)
else
    runtests(ParticleInCell2, tags = Symbol(ARGS[1]))
end
