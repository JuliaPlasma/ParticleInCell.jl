abstract type AbstractField{T, D, G} end

# TODO: Consider making these singleton types to help with dispatch
@enum GridOffset begin
    node
    edge
    face
end

struct Field{T, D, D1, G} <: AbstractField{T, D, G}
    values::Array{T, D1}
    grid::G
    offset::GridOffset

    lower_guard_cells::Int
    upper_guard_cells::Int
    index_offset::NTuple{D, Int}

    function Field(grid::G, offset::GridOffset, num_elements::Int,
        lower_guard_cells=0, upper_guard_cells=lower_guard_cells+1,
        T=Float64) where {D, G <: AbstractGrid{D}}

        index_offset = ntuple(_ -> lower_guard_cells + 1, D)
        values = zeros(T,
            (grid.num_cells .+ lower_guard_cells .+ upper_guard_cells)...,
            num_elements)

        new{T, D, D + 1, G}(values, grid, offset, lower_guard_cells,
            upper_guard_cells, index_offset)
    end
end

num_elements(f::Field) = last(size(f.values))

function Base.show(io::IO, f::Field{T, D, G}) where {T, D, G}
    print(io, "Field(grid=$(f.grid), offset=$(f.offset),
    num_elements=$(num_elements(f)), T=$T")
end

@inline function Base.eachindex(f::Field)
    return CartesianIndices(Tuple(
        UnitRange.(f.index_offset, f.index_offset .+ f.grid.num_cells .- 1)
        ))
end

@inline function cell_index_to_cell_coords(f::Field, I)
    return Tuple(I) .- f.index_offset
end

@inline function cell_coords_to_cell_index(f::Field, idxs)
    return CartesianIndex(idxs .+ f.index_offset)
end

@inline function cell_index_to_phys_coords(f::Field, I, offset=node,
    component::Int=1)
    cell_coords = cell_index_to_cell_coords(f, I)
    return cell_coords_to_phys_coords(f.grid, cell_coords, offset, component)
end

@inline function phys_coords_to_cell_index_ittr(f::Field, xs, width::Int)
    # TODO: Remove this first Tuple when Species starts using StaticArrays
    cell_coords = Tuple(phys_coords_to_cell_coords(f.grid, xs))
    cell_index = Tuple(cell_coords_to_cell_index(f, cell_coords))
    # Need to add one to lower bound to account for the fact that the particle
    # is 'in' the cell at cell_index.
    low_bounds  = cell_index .- width .+ 1
    high_bounds = cell_index .+ width
    return CartesianIndices(Tuple(UnitRange.(low_bounds, high_bounds)))
end

@inline function interp_dist(f::Field, I::CartesianIndex, xs)
    grid_pos = cell_index_to_phys_coords(f, I)
    return (xs .- grid_pos) ./ cell_lengths(f.grid)
end
