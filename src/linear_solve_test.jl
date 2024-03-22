@testitem "LinearSolve" tags = [:unit] begin
    using Statistics

    @testset "1D periodic" begin
        function test_periodic_laplacian(L, N, k, stencil, K2)
            dx = L / N
            grid = UniformCartesianGrid((0.0,), (L,), (N,), (true,))

            source_field = Field(grid, NodeOffset(), 1)
            solve_field = Field(grid, NodeOffset(), 1)

            step = LinearSolveStep(solve_field, source_field, stencil ./ ((L / N)^2))

            source_field.values .= cos.(k * range(0, L, length = N + 1))

            ParticleInCell.step!(step)

            # The matrix for linear solve has a null space, which can be seen by
            # noting that the linear equation will still be satisfied if a constant
            # value is added to each element of the solution vector (equivalent to
            # shifting the zero-point of the electric potential). In practice this
            # is not a problem because the physical quantity, the electric field,
            # is invariant in this null space. However, it does present somewhat of
            # a problem for testing the field solve in isolation. Our solution here
            # is to prescribe that the average of the potential at each grid point
            # in the domain must be zero.
            solve_field.values[end] = solve_field.values[1]
            solve_field.values .-= mean(solve_field.values[1:N])


            @test isapprox(solve_field.values, source_field.values / K2(k, dx))
        end

        L = 1.0
        for N in [16, 128]
            # kn = N/2 corresponds to the Nylquist limit
            for kn = 1:N/2
                # Three-point Laplacian
                stencil_three_point = [-1, 2, -1]

                # See Birdsall and Langdon, Eq. 2.5.13. The extra factor of pi is to
                # account for the fact that Julia computes a normalized sinc function.
                K2_three_point(k, dx) = k^2 * sinc(k * dx / 2 / pi)^2

                test_periodic_laplacian(
                    L,
                    N,
                    kn * 2pi / L,
                    stencil_three_point,
                    K2_three_point,
                )

                # Five-point Laplacian
                stencil_five_point = [-1 / 6, -1 / 3, 1, -1 / 3, -1 / 6]
                K2_five_point(k, dx) = k^2 * sinc(k * dx / 2 / pi)^2 * (2 + cos(k * dx)) / 3

                test_periodic_laplacian(
                    L,
                    N,
                    kn * 2pi / L,
                    stencil_five_point,
                    K2_five_point,
                )
            end
        end
    end

    @testset "1D grounded boundaries" begin
        grid = UniformCartesianGrid((0.0,), (1.0,), (16,), (false,))
        source_field = Field(grid, NodeOffset(), 1, 0, 0)
        solve_field = Field(grid, NodeOffset(), 1, 0, 0)

        step = LinearSolveStep(solve_field, source_field, [-1, 2, -1])

        source_field[:] .= 0.0
        source_field[end] = 5.0
        step!(step)
        @test solve_field[begin] == 0
        @test solve_field[end] == 5
    end
end
