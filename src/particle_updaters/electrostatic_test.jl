@testitem "ElectrostaticParticlePush" tags = [:unit] begin
    using StaticArrays

    g = UniformCartesianGrid((0.0,), (1.0,), (16,), (true,))
    E = Field(g, ParticleInCell2.node, 1, 1)
    E.values .= 0.1

    positions = fill(0.5, 1, 1)
    momentums = fill(0.0, 1, 1)
    weights = [1.0]
    charge = 1.0
    mass = 1.0
    species = Species(positions, momentums, weights, charge, mass)

    dt = 1.0
    step = ElectrostaticParticlePush(species, E, dt)
    step!(step)
    @test species.positions[1][1] == 0.5
    @test species.momentums[1][1] == 0.1
end
