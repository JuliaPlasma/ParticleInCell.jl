using ParticleInCell2
using Colors
using GLMakie

epsilon_0 = 8.8e-12
elec_charge = 1.6e-19
elec_mass = 9e-31

num_cells = 200
nom_density = 1e14

beam_vel = 7.78e6
particles_per_cell = 100

perturb_vel = 0.001 * beam_vel
plasma_freq = sqrt(nom_density * elec_charge^2 / elec_mass / epsilon_0)
k_int = plasma_freq / beam_vel
wavelength = 2pi / k_int

sim_length = wavelength * 3
periodic = true
grid = UniformCartesianGrid((0.0,), (sim_length,), (num_cells,), (periodic,));

field_dimension = 1
lower_guard_cells = 10
rho = Field(grid, ParticleInCell2.node, field_dimension, lower_guard_cells)
phi = Field(grid, ParticleInCell2.node, field_dimension, lower_guard_cells)
Eedge = Field(grid, ParticleInCell2.edge, field_dimension, lower_guard_cells)
Enode = Field(grid, ParticleInCell2.node, field_dimension, lower_guard_cells);

num_particles_per_stream = div(num_cells * particles_per_cell, 2)
particles_per_macro = nom_density * sim_length / num_particles_per_stream

elec_mass = 9e-31
positions = collect(0:num_particles_per_stream-1) ./ num_particles_per_stream .* sim_length
append!(positions, positions)
momentums = fill(particles_per_macro * elec_mass * beam_vel, num_particles_per_stream)
append!(momentums, -1 .* momentums)
momentums .+= (particles_per_macro * elec_mass * perturb_vel) .* sin.(positions .* k_int)

electrons = ParticleInCell2.electrons(positions, momentums, particles_per_macro);

dt = 5e-12

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

function unwrap_phase_space(species)
    points = Vector{Point2f}(undef, length(species.positions))
    for i in eachindex(species.positions)
        points[i] = Point2f(species.positions[i][1], species.momentums[i][1])
    end

    return points
end

function phase_space_colors(species)
    logocolors = Colors.JULIA_LOGO_COLORS
    colors = Vector{RGBA{Colors.N0f8}}(undef, length(species.positions))
    for i in eachindex(species.positions)
        if species.positions[i][1] < (0.5 - 1 / 6) * sim_length
            colors[i] = logocolors.red
        elseif species.positions[i][1] < (0.5 + 1 / 6) * sim_length
            colors[i] = logocolors.green
        else
            colors[i] = logocolors.purple
        end
    end

    return colors
end

points = Observable(unwrap_phase_space(electrons))
colors = Observable(phase_space_colors(electrons))

fig = Figure(resolution = (800, 400))

ylim = 3e-14
yoffset = 1e-14
ax1 = Axis(fig[1, 1], limits = (0, sim_length, -1 * ylim - yoffset, ylim))
scatter!(ax1, points, color = colors)
lbl = text!(
    ax1,
    sim_length / 2,
    -1 * ylim - yoffset,
    text = "ParticleInCell.jl",
    align = (:center, :bottom),
    font = "JuliaMono SemiBold",
    fontsize = 72,
)
hidespines!(ax1)
hidedecorations!(ax1)

record(fig, "logo.gif", 1:200, framerate = 12) do frame
    @show frame
    for i = 1:250
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

    points[] .= unwrap_phase_space(electrons)
    colors[] .= phase_space_colors(electrons)
    notify(points)
    notify(colors)
end
