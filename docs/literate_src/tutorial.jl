# # Tutorial
#md # !!! tip
#md #     You can follow along with this tutorial using the
#md #     [Jupyter notebook](tutorial.ipynb) or the plain Julia
#md #     [script](tutorial.jl).
# This tutorial will get you up and running using `ParticleInCell2` by showing
# you how to model one of the simplest phenomena in plasma physics: an
# electrostatic (or Langmuir) oscillation.
#
# A Langmuir oscillation occurs when a slab of charge in a uniform plasma is
# displaced. The resulting charge density gradient creates a restoring force
# that causes the displaced slab of charge to return to its original position.
# But---just as in a classical pendulum oscillation---the momentum of the
# charge carries it past its equilibrium point, creating an opposite charge
# gradient, and a restoring force in the opposite direction. As a result, the
# slab of charge oscillates around its equilibrium forever (at least in this
# idealized model that ignores possible damping mechanisms).
#
# To start using `ParticleInCell2`, we must first install and load the package.
# As the package is not currently registered in Julia's General registry, we
# will using `Pkg`s `develop` function. We will also need the `StaticArrays`
# package, so we load that now.
using Pkg
Pkg.develop(url = "https://github.com/adamslc/ParticleInCell2.jl")
using ParticleInCell2
using StaticArrays

# First we set some parameters of the simulation.
nom_density = 1e14
sim_length = 1.0
num_cells = 32
dx = sim_length / num_cells
num_particles = num_cells * 10
particles_per_macro = nom_density * sim_length / num_particles
dt = 1.11e-11 * 2
n_steps = 1000
k = 1
amplitude = 1e3
thermal_amp = 0.0

epsilon_0 = 8.8e-12
charge = 1.6e-19 * particles_per_macro
mass = 9e-31 * particles_per_macro

# Next, we set up the required, grid, fields, and electron species to do
# an electrostatic PIC simulation. In this example, we only consider mobile
# electrons. In order to seed a Langmuir oscillation, we give the electrons a
# sinusoidal velocity perturbation.
g = UniformCartesianGrid((0.0,), (sim_length,), (num_cells,), (true,))
rho = Field(g, ParticleInCell2.node, 1, 1)
phi = Field(g, ParticleInCell2.node, 1, 1)
Eedge = Field(g, ParticleInCell2.edge, 1, 1)
Enode = Field(g, ParticleInCell2.node, 1, 1)

positions = Vector{SVector{1,Float64}}(undef, num_particles)
momentums = Vector{SVector{1,Float64}}(undef, num_particles)
forces = Vector{SVector{1,Float64}}(undef, num_particles)
weights = Vector{Float64}(undef, num_particles)

for i in eachindex(positions)
    positions[i] = SVector((i - 1) / num_particles)
    momentums[i] = SVector(
        mass * (amplitude * sin((i - 1) / num_particles * k * 2pi) + thermal_amp * randn()),
    )
    forces[i] = SVector(0.0)
    weights[i] = 1.0
end
electrons = Species(positions, momentums, forces, weights, charge, mass)

# In the final step of the setup, we create all of the simulation steps
# required to do the electrostatic simulation.
charge_interp = BSplineChargeInterpolation(electrons, rho, 1)
comm_rho = CommunicateGuardCells(rho, true)
field_solve = PoissonSolveFFT(rho, phi)
comm_phi = CommunicateGuardCells(phi)
finite_diff = FiniteDifferenceToEdges(phi, Eedge)
comm_Eedge = CommunicateGuardCells(Eedge)
elec_ave = AverageEdgesToNodes(Eedge, Enode)
comm_Enode = CommunicateGuardCells(Enode)
elec_interp = BSplineFieldInterpolation(electrons, Enode, 1)
push = SimpleParticlePush(electrons, dt)
comm_electrons = CommunicateSpecies(electrons, g)

# Now we are ready to run the simulation. At each step, we calculate the
# current electric field energy, which will oscillate as the electrons move
# in and out of equilibrium.
electric_field_energy = Vector{Float64}(undef, n_steps)
for n = 1:n_steps
    electric_field_energy[n] = 0
    for I in eachindex(Enode)
        electric_field_energy[n] += (dx * epsilon_0 / 2) * (Enode.values[I])^2
    end

    ## TODO
    rho.values .= 0
    for i in eachindex(electrons.forces)
        electrons.forces[i] = SVector(0.0)
    end

    step!(charge_interp)
    step!(comm_rho)
    step!(field_solve)
    step!(comm_phi)
    step!(finite_diff)
    step!(comm_Eedge)
    step!(elec_ave)
    step!(comm_Enode)
    step!(elec_interp)
    step!(push)
    step!(comm_electrons)
end

# We can now use the electric field energy to calculate the plasma frequency.
using FFTW

function find_max_freq(xs, dt)
    amps = abs.(fft(xs))
    amps[1] = 0
    max_index = findmax(amps)[2]
    return fftfreq(length(xs), 1 / dt)[max_index] * 2pi
end

max_freq = abs(find_max_freq(electric_field_energy, dt))
## Divide by 2 because the electric field energy goes through a maximum twice
## per plasma oscillation
plasma_freq = max_freq / 2
@show plasma_freq

# Finally, we can compare this to the theoretically expected result:
sqrt(nom_density * charge^2 / mass / epsilon_0)
