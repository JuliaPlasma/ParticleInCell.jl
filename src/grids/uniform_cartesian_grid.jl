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

@inline function cell_coords_to_phys_coords(grid::UniformCartesianGrid, idxs)
    return grid.lower_bounds .+ idxs .* cell_lengths(grid)
end

@inline function phys_coords_to_cell_coords(grid::UniformCartesianGrid, xs)
    return (xs .- grid.lower_bounds) ./ cell_lengths(grid)
end
