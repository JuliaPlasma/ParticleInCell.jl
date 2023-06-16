@testitem "UniformCartesianGrid" begin
    @testset "cell_coords_to_phys_coords" begin
        grid = UniformCartesianGrid(
            (0.0, 0.0, 0.0),
            (1.0, 1.0, 1.0),
            (10, 10, 10),
            (false, false, false),
        )

        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (0, 0, 0)) == (0, 0, 0)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (0, 0, 0), node) == (0, 0, 0)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (0, 0, 0), edge, 1) ==
              (0.05, 0, 0)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (0, 0, 0), face, 1) ==
              (0, 0.05, 0.05)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (0, 0, 0), edge, 3) ==
              (0, 0, 0.05)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (0, 0, 0), face, 3) ==
              (0.05, 0.05, 0)

        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (5, 5, 5)) == (0.5, 0.5, 0.5)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (5, 5, 5), node) ==
              (0.5, 0.5, 0.5)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (5, 5, 5), edge, 1) ==
              (0.55, 0.5, 0.5)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (5, 5, 5), face, 1) ==
              (0.5, 0.55, 0.55)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (5, 5, 5), edge, 3) ==
              (0.5, 0.5, 0.55)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (5, 5, 5), face, 3) ==
              (0.55, 0.55, 0.5)

        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (10, 10, 10)) == (1, 1, 1)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (10, 10, 10), node) ==
              (1, 1, 1)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (10, 10, 10), edge, 1) ==
              (1.05, 1, 1)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (10, 10, 10), face, 1) ==
              (1, 1.05, 1.05)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (10, 10, 10), edge, 3) ==
              (1, 1, 1.05)
        @test ParticleInCell2.cell_coords_to_phys_coords(grid, (10, 10, 10), face, 3) ==
              (1.05, 1.05, 1)
    end

    @testset "phys_coords_to_cell_coords" begin
        grid1 = UniformCartesianGrid((0.0, 0.0), (1.0, 1.0), (10, 10), (false, false))
        @test ParticleInCell2.phys_coords_to_cell_coords(grid1, (0, 0)) == (0, 0)
        @test ParticleInCell2.phys_coords_to_cell_coords(grid1, (0.5, 0.5)) == (5, 5)
        @test ParticleInCell2.phys_coords_to_cell_coords(grid1, (1.0, 1.0)) == (10, 10)

        grid2 = UniformCartesianGrid((0.0, 0.0), (1.0, 1.0), (10, 10), (true, true))
        @test ParticleInCell2.phys_coords_to_cell_coords(grid1, (0, 0)) == (0, 0)
        @test ParticleInCell2.phys_coords_to_cell_coords(grid1, (1.0, 1.0)) == (10, 10)
        @test all(
            isapprox.(
                ParticleInCell2.phys_coords_to_cell_coords(grid1, (1.15, 1.25)),
                (11.5, 12.5),
            ),
        )
    end
end
