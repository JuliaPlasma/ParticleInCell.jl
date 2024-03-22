@testitem "BorisParticlePush" tags = [:unit] begin
    using StaticArrays
    using FFTW

    @testset "2D" begin
        g = UniformCartesianGrid((0.0,), (1.0,), (16,), (true,))
        E = Field(g, NodeOffset(), 1, 1)
        B = Field(g, NodeOffset(), 1, 1)

        positions = [SVector(0.5)]
        momentums = [SVector(0.01, 0.0)]
        weights = [1.0]
        charge = 1.0
        mass = 1.0
        species = VariableWeightSpecies(positions, momentums, weights, charge, mass)

        dt = 0.1
        push = BorisParticlePush(species, E, B, dt)

        # Check that energy is conserved with only a magnetic field
        E.values .= 0.0
        B.values .= 1.0
        initial_energy = sum(species.momentums[1] .^ 2)
        for _ = 1:10
            step!(push)
        end
        final_energy = sum(species.momentums[1] .^ 2)
        @test isapprox(initial_energy, final_energy)

        # Check that cyclotron frequency is correct
        E.values .= 0.0
        B.values .= 1.0
        species.positions .= [SVector(0.5)]
        species.momentums .= [SVector(0.01, 0.0)]

        nsteps = 200
        x = Vector{Float64}(undef, nsteps)
        for i = 1:nsteps
            x[i] = species.positions[1][1]
            step!(push)
        end

        amps = abs.(fft(x))
        amps[1] = 0
        max_index = findmax(amps)[2]
        cyclotron_freq = fftfreq(nsteps, 1 / dt)[max_index] * 2pi
        @test isapprox(cyclotron_freq, 1.0, rtol = 0.1)
    end
end
