using ParticleInCell2
using Test
using StaticArrays

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
positions = Vector{SVector{2, Float64}}(undef, n_particles)
for i in 1:n_particles
    positions[i] = SVector(rand(2)...)
end
momentums = fill(SVector(0., 0.), n_particles)
forces =    fill(SVector(0., 0.), n_particles)
weights =   fill(0., n_particles)
s = Species(positions, momentums, forces, weights, 1., 1.)

b = BSplineChargeInterpolation(s, rho, 1)
step!(1, b)

fs = PoissonSolveFFT(rho, phi)
step!(1, fs)

f = Field(g, ParticleInCell2.node, 1, 1)
comm = CommunicateGuardCells(f)
step!(1, comm)
