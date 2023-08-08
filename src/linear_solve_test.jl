@testitem "LinearSolve" tags = [:unit] begin
    @testset "1D" begin
        grid = UniformCartesianGrid((0.0,), (1.0,), (16,), (true,))
        source_field = Field(grid, ParticleInCell2.node, 1)
        solve_field = Field(grid, ParticleInCell2.node, 1)

        step = LinearSolve(nodal_field, edge_field, [-1, 2, -1])

        source_field[:] .= 0.0
        step!(step)
        @test all(solve_field[:] .== 0)
    end
end
