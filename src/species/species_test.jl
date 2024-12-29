@testitem "Species utilities" tags = [:unit] begin
    using StaticArrays
    using Distributions

    @test electrons([1.0 2.0 3.0], [4.0 5.0 6.0], [1.0, 2.0, 3.0]) ==
          VariableWeightSpecies{1,1,Float64}(
        SVector{1,Float64}.([1.0, 2.0, 3.0]),
        SVector{1,Float64}.([4.0, 5.0, 6.0]),
        [1.0, 2.0, 3.0],
        1.602176634e-19,
        9.1093837015e-31,
    )

    @test electrons([1.0 2.0 3.0], [4.0 5.0 6.0], 1.0) ==
          VariableWeightSpecies{1,1,Float64}(
        SVector{1,Float64}.([1.0, 2.0, 3.0]),
        SVector{1,Float64}.([4.0, 5.0, 6.0]),
        [1.0, 1.0, 1.0],
        1.602176634e-19,
        9.1093837015e-31,
    )

    @test electrons([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [1.0, 2.0, 3.0]) ==
          VariableWeightSpecies{1,1,Float64}(
        SVector{1,Float64}.([1.0, 2.0, 3.0]),
        SVector{1,Float64}.([4.0, 5.0, 6.0]),
        [1.0, 2.0, 3.0],
        1.602176634e-19,
        9.1093837015e-31,
    )

    @test electrons([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], 1.0) ==
          VariableWeightSpecies{1,1,Float64}(
        SVector{1,Float64}.([1.0, 2.0, 3.0]),
        SVector{1,Float64}.([4.0, 5.0, 6.0]),
        [1.0, 1.0, 1.0],
        1.602176634e-19,
        9.1093837015e-31,
    )

    dist = Normal(0.0, 1.0)
    @test ParticleInCell.quiet_velocities(4, dist) ≈
          [-1.1503493803760087, -0.3186393639643751, 0.3186393639643751, 1.1503493803760087]

    # See Birdsall and Langdon, Table 16.5a
    @test ParticleInCell.bit_reversed_sequence(6) == [0, 1 / 2, 1 / 4, 3 / 4, 1 / 8, 5 / 8]
    @test ParticleInCell.bit_reversed_sequence(6, 3) ==
          [0, 1 / 3, 2 / 3, 1 / 9, 4 / 9, 7 / 9]
end
