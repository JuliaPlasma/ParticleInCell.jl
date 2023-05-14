abstract type AbstractGrid{D} end

# A few notes on indexing conventions:
# 1. There are three different numbering systems that can refer to a location
#    in the simulation domain:
#    a. The 'phys_coords' of a point are the physical (dimensionalful)
#       coordinates associated with that point. This value can range from the
#       lower bounds to the upper bounds of the simulation. This value will
#       typically take the form Vector{T} or NTuple{N, T} where T <: Real.
#    b. The 'cell_coords' of a point is the nondimensional location of the point
#       in units of cell lengths. This value can range from 0 to num_cells - 1,
#       or outside this range if guard cells are incuded. The value will
#       typically have the type NTuple{N, Int}.
#    c. The 'cell_index' of a point is the CartesianIndex that can be used to
#       index into field arrays at that point. This value must strictly be
#       confined to axes(field.values), which, for any given dimension, will
#       typically range from 1 to num_cells + 2*num_guard_cells + 1.
# 2. The first two types of indexing, phys_coords and cell_coords, are
#    independent of the number of guard cells in a given field, and depend only
#    on grid quanties. Thus utilities for converting between these systems live
#    in this file. On the other hand, the utilities for cell_index are specific
#    the field being used, and thus those are defined where the AbstractField
#    types are defined.

struct UniformCartesianGrid{D,T,U} <: AbstractGrid{D}
    lower_bounds::NTuple{D,T}
    upper_bounds::NTuple{D,T}
    num_cells::NTuple{D,U}
    periodic::NTuple{D,Bool}
end

# Modified from WaterLily.jl
unit_vec(i, ::Val{D}) where {D} = ntuple(j -> j == i ? 1 : 0, D)
orth_vec(i, ::Val{3}) = ntuple(j -> j == i ? 0 : 1, 3)
orth_vec(::Any, ::Union{Val{1},Val{2}}) = error("orth_vec only makes sense in 3D")

@inline function cell_lengths(grid::UniformCartesianGrid)
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
        return node_coords .+ cell_lengths(grid) ./ 2 .* unit_vec(component, D)
    elseif offset == face
        return node_coords .+ cell_lengths(grid) ./ 2 .* orth_vec(component, D)
    end
end

@inline function phys_coords_to_cell_coords(grid::UniformCartesianGrid, xs)
    return floor.(Int, (xs .- grid.lower_bounds) ./ cell_lengths(grid))
end
