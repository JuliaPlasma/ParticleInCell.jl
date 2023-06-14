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

@testitem "b-spline basics" begin
    import ParticleInCell2:
        poly_integrate, binomial_expand, compute_knots, compute_bspline_coeffs

    @test poly_integrate([1.0]) == [0, 1]
    @test poly_integrate([3.0, 4.0]) == [0, 3, 2]

    @test binomial_expand([0.0, 1], 0.5) == [0.5, 1.0]
    @test binomial_expand([1.0, 0.0, 1], 0.5) == [1.25, 1.0, 1.0]

    @test compute_knots(0) == [-0.5, 0.5]
    @test compute_knots(1) == [-1, 0, 1]
    @test compute_knots(2) == [-1.5, -0.5, 0.5, 1.5]
    @test compute_knots(3) == [-2, -1, 0, 1, 2]

    # Polynomial coefficients can be computed from
    # https://en.wikipedia.org/wiki/Irwin–Hall_distribution#Special_cases
    # by translating the basis functions so that they are centered at x=0.
    @test isapprox(compute_bspline_coeffs(0), [[1]])
    @test isapprox(compute_bspline_coeffs(1), [[1, 1], [1, -1]])
    @test isapprox(
        compute_bspline_coeffs(2),
        [[9 / 8, 3 / 2, 1 / 2], [3 / 4, 0, -1], [9 / 8, -3 / 2, 1 / 2]],
    )
    @test isapprox(
        compute_bspline_coeffs(3),
        [
            [4 / 3, 2, 1, 1 / 6],
            [2 / 3, 0, -1, -1 / 2],
            [2 / 3, 0, -1, 1 / 2],
            [4 / 3, -2, 1, -1 / 6],
        ],
    )
    @test isapprox(
        compute_bspline_coeffs(4),
        [
            [625 / 384, 125 / 48, 25 / 16, 5 / 12, 1 / 24],
            [55 / 96, -5 / 24, -5 / 4, -5 / 6, -1 / 6],
            [115 / 192, 0, -5 / 8, 0, 1 / 4],
            [55 / 96, 5 / 24, -5 / 4, 5 / 6, -1 / 6],
            [625 / 384, -125 / 48, 25 / 16, -5 / 12, 1 / 24],
        ],
    )
end
