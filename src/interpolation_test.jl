@testitem "BSplineChargeInterpolation" begin
    using StaticArrays

    g = UniformCartesianGrid((0.0,), (1.0,), (10,), (true,))
    rho = Field(g, ParticleInCell2.node, 1)
    phi = Field(g, ParticleInCell2.node, 1)
    s = Species([SVector(0.5)], [SVector(0.0)], [SVector(0.0)], [1.0], 1.0, 1.0)

    bs_interp = BSplineChargeInterpolation(s, rho, 1)
    step!(bs_interp)
    @test sum(rho.values) == 1
end
# using ParticleInCell2
# using Test
# using StaticArrays

# dimension = 2
# lower_bounds = ntuple(x -> 0., dimension)
# upper_bounds = ntuple(x -> 1., dimension)
# num_cells = ntuple(x -> 16, dimension)
# periodic = ntuple(x -> true, dimension)

# g = UniformCartesianGrid(lower_bounds, upper_bounds, num_cells, periodic)
# rho = Field(g, ParticleInCell2.node, 1)
# phi = Field(g, ParticleInCell2.node, 1)


# n_particles = 100
# positions = Vector{SVector{2, Float64}}(undef, n_particles)
# for i in 1:n_particles
#     positions[i] = SVector(rand(2)...)
# end
# momentums = fill(SVector(0., 0.), n_particles)
# forces =    fill(SVector(0., 0.), n_particles)
# weights =   fill(0., n_particles)
# s = Species(positions, momentums, forces, weights, 1., 1.)

# b = BSplineChargeInterpolation(s, rho, 1)
# step!(b)

# fs = PoissonSolveFFT(rho, phi)
# step!(fs)

# f = Field(g, ParticleInCell2.node, 1, 1)
# comm = CommunicateGuardCells(f)
# step!(comm)
