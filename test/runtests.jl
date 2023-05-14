using ParticleInCell2
using Test

# @testset "ParticleInCell2.jl" begin
#     # Write your tests here.
# end

dimension = 2
lower_bounds = ntuple(x -> 0., dimension)
upper_bounds = ntuple(x -> 1., dimension)
num_cells = ntuple(x -> 16, dimension)
periodic = ntuple(x -> true, dimension)

g = UniformCartesianGrid(lower_bounds, upper_bounds, num_cells, periodic)
rho = Field(g, ParticleInCell2.node, 1)
phi = Field(g, ParticleInCell2.node, 1)


n_particles = 100
positions = rand(dimension, n_particles)
momentums = fill(0., dimension, n_particles)
forces =    fill(0., dimension, n_particles)
weights =   fill(1., dimension, n_particles)
s = Species(positions, momentums, forces, weights, 1., 1.)

b = BSplineChargeInterpolation(s, rho, 1)
step!(1, b)

fs = PoissonSolveFFT(rho, phi)
step!(1, fs)

f = Field(g, ParticleInCell2.node, 1, 1)
comm = CommunicateGuardCells(f)
step!(1, comm)
