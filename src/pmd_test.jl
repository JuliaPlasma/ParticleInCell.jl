@testitem "PMD" tags = [:unit] begin
    using Test, ParticleInCell

    grid = UniformCartesianGrid((0.0, 0.0), (1.0, 1.0), (10, 10), (true, true))

    positions = rand(2, 10)
    momentums = rand(2, 10)
    weights = rand(10)
    electrons = ParticleInCell.electrons(positions, momentums, weights)

    sim, _ = create_electrostatic_simulation(grid, [electrons], ["electrons"], 0.1)

    ParticleInCell.dump_pmd(sim, "test", 0)
end
