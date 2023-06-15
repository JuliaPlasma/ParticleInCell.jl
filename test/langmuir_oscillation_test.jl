@testitem "Langmuir Oscillation" begin
    using StaticArrays
    using FFTW

    function find_max_freq(xs, dt)
        amps = abs.(fft(xs))
        amps[1] = 0
        max_index = findmax(amps)[2]
        return fftfreq(length(xs), 1 / dt)[max_index] * 2pi
    end

    function compute_plasma_freq(nom_density)
        epsilon_0 = 8.8e-12
        charge = 1.6e-19
        mass = 9e-31

        return sqrt(nom_density * charge^2 / mass / epsilon_0)
    end

    function simulate_plamsa_freq(nom_density)
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
                mass * (
                    amplitude * sin((i - 1) / num_particles * k * 2pi) +
                    thermal_amp * randn()
                ),
            )
            forces[i] = SVector(0.0)
            weights[i] = 1.0
        end
        electrons = Species(positions, momentums, forces, weights, charge, mass)

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

        electric_field_energy = Vector{Float64}(undef, n_steps)
        for n = 1:n_steps
            electric_field_energy[n] = 0
            for I in eachindex(Enode)
                electric_field_energy[n] += (dx * epsilon_0 / 2) * (Enode.values[I])^2
            end

            # TODO
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

        max_freq = abs(find_max_freq(electric_field_energy, dt))
        # Divide by 2 because the electric field energy goes through a maximum twice
        # per plasma oscillation
        return max_freq / 2
    end

    # Allow a 10% error
    @test isapprox(compute_plasma_freq(1e14), simulate_plamsa_freq(1e14), rtol=0.1)
    @test isapprox(compute_plasma_freq(1e15), simulate_plamsa_freq(1e15), rtol=0.1)
end
