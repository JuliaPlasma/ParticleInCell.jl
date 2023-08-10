using ParticleInCell2
using BenchmarkTools
using StaticArrays

# Set up some fields
dimension = 2
lower_bounds = ntuple(x -> 0.0, dimension)
upper_bounds = ntuple(x -> 1.0, dimension)
num_cells = ntuple(x -> 16, dimension)
periodic = ntuple(x -> true, dimension)
g = UniformCartesianGrid(lower_bounds, upper_bounds, num_cells, periodic)
rho = Field(g, ParticleInCell2.node, 1)
phi = Field(g, ParticleInCell2.node, 1)
Eedge = Field(g, ParticleInCell2.edge, 2)
Enode = Field(g, ParticleInCell2.node, 2)

# Create a species
n_particles = 100
positions = Vector{SVector{2,Float64}}(undef, n_particles)
for i = 1:n_particles
    positions[i] = SVector(rand(2)...)
end
momentums = fill(SVector(0.0, 0.0), n_particles)
weights = fill(0.0, n_particles)
species = VariableWeightSpecies(positions, momentums, weights, 1.0, 1.0)

const SUITE = BenchmarkGroup()

bs_charge = BSplineChargeInterpolation(species, rho, 1)
SUITE["charge_dep"] = BenchmarkGroup(["interpolation", "particle", "field"])
SUITE["charge_dep"]["creation"] =
    @benchmarkable BSplineChargeInterpolation($species, $rho, 1)
SUITE["charge_dep"]["step"] = @benchmarkable step!($bs_charge)

fs = PoissonSolveFFT(rho, phi)
SUITE["fft_field_solve"] = BenchmarkGroup(["interpolation", "particle", "field"])
SUITE["fft_field_solve"]["creation"] = @benchmarkable PoissonSolveFFT($rho, $phi)
SUITE["fft_field_solve"]["step"] = @benchmarkable step!($fs)
