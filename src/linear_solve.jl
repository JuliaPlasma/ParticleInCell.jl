"""
Solves the linear equation `A x = b` where `x` and `b` are fields.
"""
struct LinearSolveStep{F,M,T,S} <: AbstractSimulationStep
    x::F
    b::F
    A::M
    matrix_size::T
    linear_solver::S
end

function LinearSolveStep(x::F, b::F, stencil) where {F}
    @assert x.grid === b.grid

    if !all(x.grid.periodic)
        @warn "Assuming grounded conductors for all non-periodic dimensions" grid = x.grid
    end

    # TODO: this should be part of the grid interface
    dims = x.grid.num_cells

    @assert ndims(stencil) == length(dims) "Stencil must be be the same dimension as dims"
    @assert length(stencil) % 2 == 1 "Stencil must have an odd number of entries"
    stencil_offset_index = CartesianIndex(-1 .- div.(size(stencil), 2))

    linear_indices = LinearIndices(dims)

    N = length(linear_indices)
    A = spzeros(N, N)

    for I in CartesianIndices(dims), J in eachindex(IndexCartesian(), stencil)
        # If I is a boundary cell of a non-periodic dimension, don't apply the stencil
        if any((Tuple(I) .== dims .|| Tuple(I) .== 1) .&& .!x.grid.periodic)
            A[linear_indices[I], linear_indices[I]] = 1
            continue
        end

        newI = I + J + stencil_offset_index
        roundedI = CartesianIndex(mod1.(Tuple(newI), dims))

        mat_row_index = linear_indices[I]
        mat_col_index = linear_indices[roundedI]
        A[mat_row_index, mat_col_index] += stencil[J]
    end

    linear_b = reshape(view(b.values, axes(b)...), N)
    linear_solver = init(LinearProblem(A, linear_b))

    return LinearSolveStep(x, b, A, N, linear_solver)
end

function step!(step::LinearSolveStep)
    step.linear_solver.b = reshape(view(step.b.values, axes(step.b)...), step.matrix_size)
    sol = solve!(step.linear_solver)
    xview = view(step.x.values, axes(step.x)...)
    xview .= reshape(sol.u, length.(axes(step.x)))
end
