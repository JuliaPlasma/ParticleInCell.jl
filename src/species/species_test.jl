@testitem "Species utilities" tags = [:unit] begin
    using StaticArrays

    @test electrons([1.0 2.0 3.0], [4.0 5.0 6.0], [1.0, 2.0, 3.0]) ==
          VariableWeightSpecies{1,1,Float64}(
        SVector{1,Float64}.([1.0, 2.0, 3.0]),
        SVector{1,Float64}.([4.0, 5.0, 6.0]),
        [1.0, 2.0, 3.0],
        1.602176634e-19,
        9.1093837015e-31,
    )
end
