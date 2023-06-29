struct FiniteDifferenceToEdges{F,NT} <: AbstractSimulationStep
    nodal_field::F
    edge_field::F

    edge_lengths::NT

    function FiniteDifferenceToEdges(nodal_field::F, edge_field::F) where {F}
        @assert nodal_field.offset == node
        @assert edge_field.offset == edge
        @assert edge_field.upper_guard_cells >= 1

        edge_lengths = cell_lengths(nodal_field.grid)

        new{F,typeof(edge_lengths)}(nodal_field, edge_field, edge_lengths)
    end
end

# TODO: document that this calculates -1*grad(nodal_field). Maybe make it an option?
function step!(step::FiniteDifferenceToEdges{F}) where {T,D,F<:AbstractField{T,D}}
    nodal_field = step.nodal_field
    edge_field = step.edge_field

    # TODO: flip the order of these loops?
    for I in eachindex(edge_field)
        for d = 1:D
            dI = CartesianIndex(unit_vec(d, Val{D}()))
            # TODO: need to update the indexing into edge field, because it has
            # multiple components
            edge_field[I] = -1 * (nodal_field[I+dI] - nodal_field[I]) / step.edge_lengths[d]
        end
    end
end

struct AverageEdgesToNodes{F} <: AbstractSimulationStep
    edge_field::F
    nodal_field::F

    function AverageEdgesToNodes(edge_field::F, nodal_field::F) where {F}
        @assert edge_field.offset == edge
        @assert nodal_field.offset == node
        @assert edge_field.lower_guard_cells >= 1
        @assert edge_field.upper_guard_cells >= 1

        new{F}(edge_field, nodal_field)
    end
end

function step!(step::AverageEdgesToNodes{F}) where {T,D,F<:AbstractField{T,D}}
    edge_field = step.edge_field
    nodal_field = step.nodal_field

    for I in eachindex(edge_field)
        for d = 1:D
            dI = CartesianIndex(unit_vec(d, Val{D}()))
            nodal_field[I] = (edge_field[I-dI] + edge_field[I]) / 2
        end
    end
end

struct CommunicateGuardCells{F} <: AbstractSimulationStep
    field::F
    sum_guards::Bool
end

CommunicateGuardCells(f) = CommunicateGuardCells(f, false)

function step!(step::CommunicateGuardCells{F}) where {T,D,F<:AbstractField{T,D}}
    field = step.field
    grid = field.grid

    # TODO: need to add an extra Colon for the final dimension of the values array
    for d = 1:D
        lower_guard_indices =
            ntuple(j -> j == d ? UnitRange(1, field.lower_guard_cells) : Colon(), D)
        lower_match_indices = ntuple(
            j ->
                j == d ?
                UnitRange(
                    grid.num_cells[d] + 1,
                    grid.num_cells[d] + field.lower_guard_cells,
                ) : Colon(),
            D,
        )
        if step.sum_guards
            field.values[lower_match_indices...] .+= field.values[lower_guard_indices...]
        end
        field.values[lower_guard_indices...] .= field.values[lower_match_indices...]

        upper_guard_indices = ntuple(
            j ->
                j == d ?
                UnitRange(
                    field.lower_guard_cells + grid.num_cells[d] + 1,
                    field.lower_guard_cells + grid.num_cells[d] + field.upper_guard_cells,
                ) : Colon(),
            D,
        )
        upper_match_indices = ntuple(
            j ->
                j == d ?
                UnitRange(
                    field.lower_guard_cells + 1,
                    field.lower_guard_cells + field.upper_guard_cells,
                ) : Colon(),
            D,
        )

        if step.sum_guards
            field.values[upper_match_indices...] .+= field.values[upper_guard_indices...]
        end
        field.values[upper_guard_indices...] .= field.values[upper_match_indices...]
    end
end