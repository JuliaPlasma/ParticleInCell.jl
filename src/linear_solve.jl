"""
Solves the linear equation `A x = b` where `x` and `b` are fields.
"""
struct LinearSolve{F,M,T,S} <: AbstractSimulationStep
    x::F
    b::F
    A::M
    matrix_size::T
    linear_solver::S
end

function LinearSolve(x::F, b::F, stencil) where {F}
    @assert x.grid === b.grid

    @assert ndims(stencil) == length(dims) "Stencil must be be the same dimension as dims"
    @assert length(stencil) % 2 == 1 "Stencil must have an odd number of entries"
    stencil_offset_index = CartesianIndex(-1 .- div.(size(stencil), 2))

    linear_indices = LinearIndices(dims)

    N = length(linear_indices)
    A = spzeros(N, N)

    for I in CartesianIndices(dims), J in eachindex(IndexCartesian(), stencil)
        newI = I + J + stencil_offset_index
        roundedI = CartesianIndex(mod1.(Tuple(newI), dims))

        mat_row_index = linear_indices[I]
        mat_col_index = linear_indices[roundedI]
        A[mat_row_index, mat_col_index] += stencil[J]
    end

    linear_b = reshape(b.values, N)
    linear_solver = init(LinearProblem(A, linear_b))

    return LinearSolve(x, b, A, N, linear_solver)
end

function step!(step::LinearSolve)
    step.linear_solver.b = reshape(step.b.values, step.N)
    sol = solve!(step.linear_solver)
    x.values .= reshape(sol.u, size(x.values))
end
