@testitem "BSplineChargeInterpolation" begin
    using StaticArrays

    function sum_charge(rho)
        cell_volume = prod(ParticleInCell2.cell_lengths(rho.grid))
        charge = 0.0
        for I in eachindex(rho)
            charge += cell_volume * rho.values[I]
        end
        return charge
    end

    g = UniformCartesianGrid((0.0,), (1.0,), (10,), (true,))
    rho = Field(g, ParticleInCell2.node, 1)
    phi = Field(g, ParticleInCell2.node, 1)
    s = Species([SVector(0.5)], [SVector(0.0)], [SVector(0.0)], [1.0], 1.0, 1.0)

    bs_interp = BSplineChargeInterpolation(s, rho, 1)
    step!(bs_interp)
    @test sum_charge(rho) == 1
end
