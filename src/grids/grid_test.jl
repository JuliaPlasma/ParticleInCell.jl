@testitem "UniformCartesianGrid" tags = [:unit] begin
    @testset "cell_coords_to_phys_coords" begin
        grid = UniformCartesianGrid(
            (0.0, 0.0, 0.0),
            (1.0, 1.0, 1.0),
            (10, 10, 10),
            (false, false, false),
        )

        @test ParticleInCell.cell_coords_to_phys_coords(grid, (0, 0, 0)) == (0, 0, 0)

        @test ParticleInCell.cell_coords_to_phys_coords(grid, (5, 5, 5)) == (0.5, 0.5, 0.5)

        @test ParticleInCell.cell_coords_to_phys_coords(grid, (10, 10, 10)) == (1, 1, 1)
    end

    @testset "phys_coords_to_cell_coords" begin
        grid1 = UniformCartesianGrid((0.0, 0.0), (1.0, 1.0), (10, 10), (false, false))
        @test ParticleInCell.phys_coords_to_cell_coords(grid1, (0, 0)) == (0, 0)
        @test ParticleInCell.phys_coords_to_cell_coords(grid1, (0.5, 0.5)) == (5, 5)
        @test ParticleInCell.phys_coords_to_cell_coords(grid1, (1.0, 1.0)) == (10, 10)

        grid2 = UniformCartesianGrid((0.0, 0.0), (1.0, 1.0), (10, 10), (true, true))
        @test ParticleInCell.phys_coords_to_cell_coords(grid1, (0, 0)) == (0, 0)
        @test ParticleInCell.phys_coords_to_cell_coords(grid1, (1.0, 1.0)) == (10, 10)
        @test all(
            isapprox.(
                ParticleInCell.phys_coords_to_cell_coords(grid1, (1.15, 1.25)),
                (11.5, 12.5),
            ),
        )
    end
end
