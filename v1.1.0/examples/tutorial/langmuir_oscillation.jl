using ParticleInCell
using CairoMakie
CairoMakie.activate!(type = "svg") #hide
docs_theme = Theme(
    fontsize = 20,
    fonts = (; regular = "Lato", bold = "Lato Bold"),
    palette = (color = [:black, :red, :blue],),
    Axis = (xgridvisible = false, ygridvisible = false),
    Scatter = (; cycle = [:color]),
    Lines = (; cycle = [:color]),
)
set_theme!(docs_theme) #hide

sim_length = 1.0
number_density = 1e14
num_macroparticles = 320
particles_per_macro = number_density * sim_length / num_macroparticles

positions = collect(0:num_macroparticles-1) ./ num_macroparticles;

k = 1 * 2pi / sim_length
amplitude = 1e3
elec_mass = 9e-31
momentums = (particles_per_macro * elec_mass * amplitude) .* sin.(positions .* k);

fig = Figure(size = (1000, 400))
ax = Axis(
    fig[1, 1],
    title = "Electron phase space",
    xlabel = "Position (m)",
    ylabel = "Momentum (arb. units)",
    limits = ((0, 1), nothing),
)
scatter!(ax, positions, momentums ./ maximum(momentums))
fig

electrons = ParticleInCell.electrons(positions, momentums, particles_per_macro);

num_cells = 32
dx = sim_length / num_cells
periodic = true
grid = UniformCartesianGrid((0.0,), (sim_length,), (num_cells,), (periodic,));

epsilon_0 = 8.8e-12
elec_charge = 1.6e-19
elec_mass = 9e-31
expected_plasma_freq = sqrt(number_density * elec_charge^2 / elec_mass / epsilon_0)
expected_plasma_period = 2pi / expected_plasma_freq

dt = 0.01 / expected_plasma_freq

sim, fields = create_electrostatic_simulation(grid, [electrons], dt);

Enode = fields[:Enode];

n_steps = 10000

electric_field_energy = Vector{Float64}(undef, n_steps)

for n = 1:n_steps
    # Calculate the electric field energy
    electric_field_energy[n] = 0
    for I in eachindex(Enode)
        electric_field_energy[n] += (dx * epsilon_0 / 2) * (Enode.values[I])^2
    end

    step!(sim)
end

normalized_times = collect(range(1, n_steps)) .* dt .* expected_plasma_freq
fig = Figure(size = (1000, 400))
ax = Axis(
    fig[1, 1],
    title = "Electric Field Energy",
    xlabel = "Normalized Time",
    ylabel = "Energy",
    limits = ((0, 100), nothing),
)
lines!(normalized_times, electric_field_energy)
fig

using FFTW

freqs = fftfreq(n_steps, 1 / dt) .* 2pi
freq_amps = abs.(fft(electric_field_energy))

fig = Figure(size = (1000, 400))
ax = Axis(
    fig[1, 1],
    title = "Electric Field Energy Spectrum",
    xlabel = "Frequency",
    ylabel = "Amplitude",
)
lines!(ax, freqs, freq_amps)
fig

cutoff_index = round(Int, n_steps * 0.005)
fig = Figure(size = (1000, 400))
ax = Axis(
    fig[1, 1],
    title = "Electric Field Energy Spectrum",
    xlabel = "Frequency",
    ylabel = "Amplitude",
)
lines!(ax, freqs[1:cutoff_index], freq_amps[1:cutoff_index])
fig

freq_amps .= ifelse.(freqs .< 5e8, 0, freq_amps)
max_index = findmax(freq_amps)[2]
max_freq = freqs[max_index]

# Divide by 2 because the electric field energy goes through a maximum twice
# per plasma oscillation, and take the absolute value because we don't care
# about the phase of the oscillation.
measured_plasma_freq = abs(max_freq / 2)

fig = Figure(size = (1000, 400))
ax = Axis(
    fig[1, 1],
    title = "Electric Field Energy Spectrum",
    xlabel = "Frequency",
    ylabel = "Amplitude",
)
lines!(ax, freqs[1:cutoff_index], freq_amps[1:cutoff_index])
vlines!(ax, [measured_plasma_freq * 2])
fig

relative_error = (measured_plasma_freq - expected_plasma_freq) / expected_plasma_freq

function measure_plasma_frequency(number_density, wavenumber, normalized_timestep = 0.01)
    sim_length = 1.0
    num_cells = 32
    dx = sim_length / num_cells

    num_macroparticles = 10 * num_cells
    particles_per_macro = number_density * sim_length / num_macroparticles

    perturb_amplitude = 1e3
    elec_mass = 9e-31
    positions = collect(0:num_macroparticles-1) ./ num_macroparticles
    momentums =
        (particles_per_macro * elec_mass * perturb_amplitude) .*
        sin.(positions .* wavenumber)

    electrons = ParticleInCell.electrons(positions, momentums, particles_per_macro)

    grid = UniformCartesianGrid((0.0,), (sim_length,), (num_cells,), (true,))

    epsilon_0 = 8.8e-12
    elec_charge = 1.6e-19
    expected_plasma_freq = sqrt(number_density * elec_charge^2 / elec_mass / epsilon_0)
    dt = normalized_timestep / expected_plasma_freq
    sim, fields = create_electrostatic_simulation(grid, [electrons], dt)
    Enode = fields[:Enode]

    # We want to simulate the same length of physical time for each simulation,
    # so we need to have more steps when the timestep is shorter.
    n_steps = round(Int, 100 / normalized_timestep)
    electric_field_energy = Vector{Float64}(undef, n_steps)

    epsilon_0 = 8.8e-12
    for n = 1:n_steps
        # Calculate the electric field energy
        electric_field_energy[n] = 0
        for I in eachindex(Enode)
            electric_field_energy[n] += (dx * epsilon_0 / 2) * (Enode.values[I])^2
        end

        step!(sim)
    end

    freqs = fftfreq(n_steps, 1 / dt) .* 2pi
    freq_amps = abs.(fft(electric_field_energy))

    freq_amps[1] = 0
    max_index = findmax(freq_amps)[2]
    max_freq = freqs[max_index]
    plasma_freq = abs(max_freq / 2)

    return plasma_freq
end;

wavenumbers = 2 * pi .* [2, 4, 8, 16]
plasma_freqs = measure_plasma_frequency.(1e14, wavenumbers, 0.001)

dx = 1 / 32
normalized_wavenumbers = wavenumbers .* dx
relative_errors = abs.(plasma_freqs ./ expected_plasma_freq .- 1)

fig = Figure(size = (1000, 400))
ax = Axis(
    fig[1, 1],
    title = "Frequency error",
    xlabel = L"k \Delta x",
    ylabel = L"\left|\frac{\omega}{\omega_p} - 1\right|",
    xscale = log10,
    yscale = log10,
)
scatter!(ax, normalized_wavenumbers, relative_errors)
fig

using LsqFit

model(x, p) = p[1] .+ x .* p[2]
fit = curve_fit(model, log10.(normalized_wavenumbers), log10.(relative_errors), [0.0, 2.0])
params = coef(fit)

lines!(
    ax,
    normalized_wavenumbers,
    (10^params[1]) .* normalized_wavenumbers .^ params[2],
    color = :red,
)
fig

corrected_plasma_freqs =
    plasma_freqs ./ (1 .- (10^params[1]) .* normalized_wavenumbers .^ params[2])
corrected_relative_errors =
    (corrected_plasma_freqs .- expected_plasma_freq) ./ expected_plasma_freq

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
