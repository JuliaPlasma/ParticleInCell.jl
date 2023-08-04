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
