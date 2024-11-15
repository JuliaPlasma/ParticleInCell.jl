struct Simulation{T}
    grid::Any
    fields::Vector
    field_names::Vector{String}
    species::Vector
    species_names::Vector{String}
    steps::Vector{AbstractSimulationStep}

    time::T
    dt::T
end

function step!(sim::Simulation)
    for step in sim.steps
        step!(step)
    end

    sim.time += sim.dt
end

function create_electrostatic_simulation(
    grid,
    species::Vector{S},
    species_names::Vector{String},
    dt;
    interpolation_order::Int = 1,
) where {D,V,S<:AbstractSpecies{D,V}}

    # Create fields
    lower_guard_cells = div(interpolation_order, 2) + 1
    rho = Field(grid, NodeOffset(), 1, lower_guard_cells)
    phi = Field(grid, NodeOffset(), 1, lower_guard_cells)
    Eedge = Field(grid, EdgeOffset(), V, lower_guard_cells)
    Enode = Field(grid, NodeOffset(), V, lower_guard_cells)

    fields = [rho, phi, Eedge, Enode]
    field_names = ["rho", "phi", "Eedge", "Enode"]

    sim = Simulation(
        grid,
        fields,
        field_names,
        species,
        species_names,
        AbstractSimulationStep[],
        zero(dt),
        dt,
    )

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

function dump_pmd(
    sim::Simulation,
    basename,
    dump_number;
    author = nothing,
    machine = gethostname(),
    comment = nothing,
)
    filename = "$(basename)_$(dump_number).h5"
    h5open(filename, "w") do f

        attributes(f)["openPMD"] = "1.1.0"
        attributes(f)["openPMDextension"] = UInt32(0)
        attributes(f)["basePath"] = "/data/%T"

        attributes(f)["meshesPath"] = "meshes/"
        attributes(f)["particlesPath"] = "particles/"

        if author !== nothing
            attributes(f)["author"] = author
        end
        attributes(f)["software"] = "ParticleInCell.jl"
        attributes(f)["softwareVersion"] = string(pkgversion(ParticleInCell))
        # This should include a timezone for full PMD compliance
        attributes(f)["date"] = string(now())

        # TODO: add package versions, and kwarg to specify more dependencies
        attributes(f)["softwareDependencies"] = "julia@$(VERSION)"
        attributes(f)["machine"] = machine

        if comment !== nothing
            attributes(f)["comment"] = comment
        end

        attributes(f)["iterationEncoding"] = "fileBased"
        attributes(f)["iterationFormat"] = "$(basename)_%T.h5"

        data_group = create_group(f, "data")
        dump_group = create_group(data_group, string(dump_number))

        attributes(dump_group)["time"] = sim.time
        attributes(dump_group)["dt"] = sim.dt
        attributes(dump_group)["timeUnitSI"] = 1.0

        meshes_group = create_group(dump_group, "meshes")

        for (field, field_name) in zip(sim.fields, sim.field_names)
            dump_pmd_mesh(field, field_name, meshes_group)
        end

        particles_group = create_group(dump_group, "particles")

        for (species, species_name) in zip(sim.species, sim.species_names)
            dump_pmd_particles(species, species_name, sim.grid, particles_group)
        end
    end

end

function dump_pmd_mesh(field, field_name, meshes_group)
    field_dimension = last(last(axes(field.values)))
    if field_dimension == 1
        write(meshes_group, field_name, view(field.values, eachindex(field)))
        mesh_item = meshes_group[field_name]
        attributes(mesh_item)["position"] = collect(position_offset(field, 1))
        attributes(mesh_item)["unitSI"] = 1.0
    else
        mesh_item = create_group(meshes_group, field_name)
        for dim = 1:field_dimension
            dim_name = axis_label(field.grid, dim)
            write(mesh_item, dim_name, view(field.values, eachindex(field), dim))
            attributes(mesh_item[dim_name])["position"] =
                collect(position_offset(field, dim))
            attributes(mesh_item[dim_name])["unitSI"] = 1.0
        end
    end

    attributes(mesh_item)["geometry"] = geometry_name(field.grid)
    attributes(mesh_item)["dataOrder"] = "F"
    attributes(mesh_item)["axisLabels"] =
        collect([axis_label(field.grid, dim) for dim = 1:field_dimension])
    # TODO: this doesn't support variable grids, but I don't think the current openPMD
    # standard has any way to support that...
    attributes(mesh_item)["gridSpacing"] = collect(cell_lengths(field.grid))
    # TODO: this only works for UniformCartesianGrid
    attributes(mesh_item)["gridGlobalOffset"] = collect(field.grid.lower_bounds)
    attributes(mesh_item)["gridUnitSI"] = 1.0

    # TODO: Need some way to store this with each field. Might make sense to incorporate
    # Unitful to make specifying units easier.
    attributes(mesh_item)["unitDimension"] = zeros(7)
    # TODO: Also need to define this in the field
    attributes(mesh_item)["timeOffset"] = 0.0
end

function dump_pmd_particles(
    species::S,
    species_name,
    grid,
    particles_group,
) where {D,V,S<:AbstractSpecies{D,V}}
    species_group = create_group(particles_group, species_name)

    position_group = create_group(species_group, "position")
    position_offset_group = create_group(species_group, "positionOffset")
    attributes(position_group)["unitDimension"] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    attributes(position_offset_group)["unitDimension"] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    attributes(position_offset_group)["timeOffset"] = 0.0

    for i = 1:D
        label = axis_label(grid, i)
        write(
            position_group,
            label,
            [particle_position(species, p)[i] for p in eachindex(species)],
        )
        attributes(position_offset_group)[label] = 0.0
        attributes(position_group[label])["unitSI"] = 1.0
    end

    momentum_group = create_group(species_group, "momentum")
    attributes(momentum_group)["unitDimension"] = [1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    attributes(momentum_group)["timeOffset"] = 0.5
    for i = 1:V
        label = axis_label(grid, i)
        write(
            momentum_group,
            label,
            [particle_momentum(species, p)[i] for p in eachindex(species)],
        )
        attributes(momentum_group[label])["unitSI"] = 1.0
    end

    write(
        species_group,
        "weight",
        [particle_weight(species, p) for p in eachindex(species)],
    )
    attributes(species_group["weight"])["unitDimension"] = zeros(7)
    attributes(species_group["weight"])["timeOffset"] = 0.0
    attributes(species_group["weight"])["unitSI"] = 1.0

    charge_group = create_group(species_group, "charge")
    attributes(charge_group)["value"] = physical_charge(species, 1)
    attributes(charge_group)["shape"] = [length(eachindex(species))]
    attributes(charge_group)["unitDimension"] = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    attributes(charge_group)["timeOffset"] = 0.0
    attributes(charge_group)["unitSI"] = 1.0

    mass_group = create_group(species_group, "mass")
    attributes(mass_group)["value"] = physical_mass(species, 1)
    attributes(mass_group)["shape"] = [length(eachindex(species))]
    attributes(mass_group)["unitDimension"] = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    attributes(mass_group)["timeOffset"] = 0.0
    attributes(mass_group)["unitSI"] = 1.0

    # TODO: implement particle patches
end
