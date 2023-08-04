@testitem "FiniteDifferenceToEdges" tags = [:unit] begin
    @testset "1D" begin
        grid = UniformCartesianGrid((0.0,), (1.0,), (10,), (false,))
        nodal_field = Field(grid, ParticleInCell2.node, 1)
        edge_field = Field(grid, ParticleInCell2.edge, 1)

        step = FiniteDifferenceToEdges(nodal_field, edge_field)
        @test step.edge_lengths == (0.1,)

        nodal_field[:] .= 1.0
        step!(step)
        @test all(edge_field[:] .== 0)

        nodal_field[:] = collect(range(0, 10))
        step!(step)
        @test edge_field[5] .== -10
        @test edge_field[11] .== 0
    end
end

@testitem "AverageEdgesToNodes" tags = [:unit] begin
    @testset "1D" begin
        grid = UniformCartesianGrid((0.0,), (1.0,), (10,), (false,))
        nodal_field = Field(grid, ParticleInCell2.node, 1, 1)
        edge_field = Field(grid, ParticleInCell2.edge, 1, 1)

        step = AverageEdgesToNodes(edge_field, nodal_field)

        edge_field.values[:] .= 2
        step!(step)
        @test nodal_field[5] == 2
    end
end
