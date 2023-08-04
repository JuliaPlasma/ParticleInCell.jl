@testitem "BorisParticlePush" tags = [:unit] begin
    using StaticArrays

    g = UniformCartesianGrid((0.0, 0.0), (1.0, 1.0), (16, 16), (true, true))
    E = Field(g, ParticleInCell2.node, 2, 1)
    B = Field(g, ParticleInCell2.node, 2, 1)

    positions = [SVector(0.5, 0.5)]
    momentums = [SVector(0.0, 1.0)]
    weights = [1.0]
    charge = 1.0
    mass = 1.0
    species = Species(positions, momentums, weights, charge, mass)

    dt = 1.0
    push = BorisParticlePush(species, E, B, dt)
end
