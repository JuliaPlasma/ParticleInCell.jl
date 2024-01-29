struct Simulation
    steps::Vector{AbstractSimulationStep}
end

function step!(sim::Simulation)
    for step in sim.steps
        step!(step)
    end
end

function create_electrostatic_simulation(grid, species::Vector{S}, dt; interpolation_order::Int=1) where {D,V,S<:AbstractSpecies{D,V}}
    sim = Simulation(AbstractSimulationStep[])

    # Create fields
    lower_guard_cells = div(interpolation_order, 2) + 1
    rho = Field(grid, ParticleInCell.node,   V, lower_guard_cells)
    phi = Field(grid, ParticleInCell.node,   V, lower_guard_cells)
    Eedge = Field(grid, ParticleInCell.edge, V, lower_guard_cells)
    Enode = Field(grid, ParticleInCell.node, V, lower_guard_cells)

    # Zero out charge density then deposit and communicate rho
    push!(sim.steps, ZeroField(rho))
    for s in species
        push!(sim.steps, BSplineChargeInterpolation(s, rho, interpolation_order))
    end
    push!(sim.steps, CommunicateGuardCells(rho, true))

    # Field solve and communicate phi
    push!(sim.steps, PoissonSolveFFT(rho, phi))
    push!(sim.steps, CommunicateGuardCells(phi))

    # Calculate Eedge and communicate
    push!(sim.steps, FiniteDifferenceToEdges(phi, Eedge))
    push!(sim.steps, CommunicateGuardCells(Eedge))

    # Calculate Enode and communicate
    push!(sim.steps, AverageEdgesToNodes(Eedge, Enode))
    push!(sim.steps, CommunicateGuardCells(Enode))

    # Push particles for each species
    for s in species
        push!(sim.steps, ElectrostaticParticlePush(s, Enode, dt))
        push!(sim.steps, CommunicateSpecies(s, grid))
    end

    return sim, (; rho, phi, Eedge, Enode)
end
