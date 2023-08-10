@testitem "VariableWeightSpecies" tags = [:unit] begin
    using StaticArrays

    species =
        VariableWeightSpecies([1.0 2.0 3.0], [4.0 5.0 6.0], [1.0, 2.0, 3.0], 5.0, 10.0)
    @test species == VariableWeightSpecies{1,1,Float64}(
        SVector{1,Float64}.([1.0, 2.0, 3.0]),
        SVector{1,Float64}.([4.0, 5.0, 6.0]),
        [1.0, 2.0, 3.0],
        5.0,
        10.0,
    )

    @test particle_charge(species, 1) == 5.0
    @test particle_charge(species, 2) == 10.0
    @test physical_charge(species, 1) == 5.0
    @test physical_charge(species, 2) == 5.0

    @test particle_mass(species, 1) == 10.0
    @test particle_mass(species, 2) == 20.0
    @test physical_mass(species, 1) == 10.0
    @test physical_mass(species, 2) == 10.0

    @test particle_weight(species, 1) == 1.0
    @test particle_weight(species, 2) == 2.0

    @test particle_position(species, 1)[1] == 1.0
    particle_position!(species, 1, SVector(2.0))
    @test particle_position(species, 1)[1] == 2.0

    @test particle_momentum(species, 1)[1] == 4.0
    particle_momentum!(species, 1, SVector(2.0))
    @test particle_momentum(species, 1)[1] == 2.0

    @test particle_velocity(species, 2)[1] == 5.0 / (10.0 * 2.0)

    @test physical_momentum(species, 2)[1] == 5.0 / 2.0
    @test physical_momentum(species, 3)[1] == 6.0 / 3.0

    @test eachindex(species) == Base.OneTo(3)
end
