@testitem "CommunicateGuardCells" tags = [:unit] begin
    @testset "1D" begin
        grid = UniformCartesianGrid((0.0,), (1.0,), (10,), (false,))
        field = Field(grid, ParticleInCell2.node, 1)

        step = CommunicateGuardCells(field)

        field[1] = 1.0
        field[11] = 0.0
        @test field[1] == 1.0
        @test field[11] == 0.0
        step!(step)
        @test field[1] == 1.0
        @test field[11] == 1.0
    end
end
