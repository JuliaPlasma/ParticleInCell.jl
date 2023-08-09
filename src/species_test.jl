@testitem "Species" tags = [:unit] begin
    using StaticArrays

    @test Species([1.0 2.0 3.0], [4.0 5.0 6.0], [1.0, 2.0, 3.0], 1.0, 1.0) ==
          Species{1,1,Float64}(
        SVector{1,Float64}.([1.0, 2.0, 3.0]),
        SVector{1,Float64}.([4.0, 5.0, 6.0]),
        [1.0, 2.0, 3.0],
        1.0,
        1.0,
    )
end
