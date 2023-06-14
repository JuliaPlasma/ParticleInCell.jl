"""
    UniformCartesianGrid(lower_bounds, upper_bounds, num_cells, periodic)

The simplest grid type, which represents a set of equally spaced rectangular cells.
The grid can optionally be periodic in one or more dimensions.
"""
struct UniformCartesianGrid{D,T,U} <: AbstractGrid{D}
    lower_bounds::NTuple{D,T}
    upper_bounds::NTuple{D,T}
    num_cells::NTuple{D,U}
    periodic::NTuple{D,Bool}
end

@inline function sim_lengths(grid::UniformCartesianGrid)
    return (grid.upper_bounds .- grid.lower_bounds)
end

@inline function cell_lengths(grid::UniformCartesianGrid, idxs = nothing)
    return (grid.upper_bounds .- grid.lower_bounds) ./ grid.num_cells
end

# TODO: inline this function once the offsets become singltons, and this
# function becomes performant
function cell_coords_to_phys_coords(
    grid::UniformCartesianGrid{D},
    idxs,
    offset = node,
    component::Int = 1,
) where {D}
    node_coords = grid.lower_bounds .+ idxs .* cell_lengths(grid)

    if offset == node
        return node_coords
    elseif offset == edge
        return node_coords .+ cell_lengths(grid) ./ 2 .* unit_vec(component, Val(D))
    elseif offset == face
        return node_coords .+ cell_lengths(grid) ./ 2 .* orth_vec(component, Val(D))
    end
end

@inline function phys_coords_to_cell_coords(grid::UniformCartesianGrid, xs)
    return floor.(Int, (xs .- grid.lower_bounds) ./ cell_lengths(grid))
end
