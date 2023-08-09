using Pkg
Pkg.develop(url = "https://github.com/adamslc/ParticleInCell2.jl")
using ParticleInCell2
using StaticArrays

sim_length = 1.0
num_cells = 32
periodic = true
grid = UniformCartesianGrid((0.0,), (sim_length,), (num_cells,), (periodic,));

field_dimension = 1
lower_guard_cells = 1
rho = Field(grid, ParticleInCell2.node, field_dimension, lower_guard_cells)
phi = Field(grid, ParticleInCell2.node, field_dimension, lower_guard_cells)
Eedge = Field(grid, ParticleInCell2.edge, field_dimension, lower_guard_cells)
Enode = Field(grid, ParticleInCell2.node, field_dimension, lower_guard_cells);

nom_density = 1e14
dx = sim_length / num_cells
num_particles = num_cells * 10
particles_per_macro = nom_density * sim_length / num_particles
k = 1
amplitude = 1e3
thermal_amp = 0.0

epsilon_0 = 8.8e-12
charge = 1.6e-19 * particles_per_macro
mass = 9e-31 * particles_per_macro

positions = Vector{SVector{1,Float64}}(undef, num_particles)
momentums = Vector{SVector{1,Float64}}(undef, num_particles)
weights = Vector{Float64}(undef, num_particles)

for i in eachindex(positions)
    positions[i] = SVector((i - 1) / num_particles)
    momentums[i] = SVector(
        mass * (amplitude * sin((i - 1) / num_particles * k * 2pi) + thermal_amp * randn()),
    )
    weights[i] = 1.0
end
electrons = Species(positions, momentums, weights, charge, mass);

dt = 1.11e-11 * 2

charge_interp = BSplineChargeInterpolation(electrons, rho, 1)
comm_rho = CommunicateGuardCells(rho, true)
field_solve = PoissonSolveFFT(rho, phi)
comm_phi = CommunicateGuardCells(phi)
finite_diff = FiniteDifferenceToEdges(phi, Eedge)
comm_Eedge = CommunicateGuardCells(Eedge)
elec_ave = AverageEdgesToNodes(Eedge, Enode)
comm_Enode = CommunicateGuardCells(Enode)
push = ElectrostaticParticlePush(electrons, Enode, dt)
comm_electrons = CommunicateSpecies(electrons, grid);

n_steps = 1000

electric_field_energy = Vector{Float64}(undef, n_steps)

for n = 1:n_steps
    # Calculate the electric field energy
    electric_field_energy[n] = 0
    for I in eachindex(Enode)
        electric_field_energy[n] += (dx * epsilon_0 / 2) * (Enode.values[I])^2
    end

    # TODO
    rho.values .= 0

    step!(charge_interp)
    step!(comm_rho)
    step!(field_solve)
    step!(comm_phi)
    step!(finite_diff)
    step!(comm_Eedge)
    step!(elec_ave)
    step!(comm_Enode)
    step!(push)
    step!(comm_electrons)
end

using CairoMakie
CairoMakie.activate!(type = "svg")

times = collect(range(1, n_steps)) .* dt
lines(
    times,
    electric_field_energy;
    axis = (; title = "Electric Field Energy", xlabel = "Time (s)", ylabel = "Energy"),
)

using FFTW

freqs = fftfreq(n_steps, 1 / dt) .* 2pi
freq_amps = abs.(fft(electric_field_energy))

lines(
    freqs,
    freq_amps;
    axis = (;
        title = "Electric Field Energy Frequency Spectrum",
        xlabel = "Frequency (1/s)",
        ylabel = "Amplitude",
    ),
)

cutoff_index = round(Int, n_steps * 0.05)
lines(
    freqs[1:cutoff_index],
    freq_amps[1:cutoff_index];
    axis = (;
        title = "Electric Field Energy Frequency Spectrum",
        xlabel = "Frequency (1/s)",
        ylabel = "Amplitude",
    ),
)

freq_amps[1] = 0
max_index = findmax(freq_amps)[2]
max_freq = freqs[max_index]

# Divide by 2 because the electric field energy goes through a maximum twice
# per plasma oscillation
plasma_freq = max_freq / 2

elec_charge = 1.6e-19
elec_mass = 9e-31
theory_plasma_freq = sqrt(nom_density * elec_charge^2 / elec_mass / epsilon_0)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

