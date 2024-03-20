abstract type AbstractField{T,N,D,G} end

# TODO: Consider making these singleton types to help with dispatch
@enum GridOffset begin
    node
    edge
    face
end

"""
    Field(grid, offset, num_elements,[ lower_guard_cells = 0,
        [upper_guard_cells = lower_guard_cells + 1]])

Stores the values of a field.
"""
struct Field{T,N,D,D1,G} <: AbstractField{T,N,D,G}
    values::Array{T,D1}
    grid::G
    offset::GridOffset

    lower_guard_cells::Int
    upper_guard_cells::Int
    index_offset::NTuple{D,Int}

    function Field(
        grid::G,
        offset::GridOffset,
        num_elements::Int,
        lower_guard_cells = 0,
        upper_guard_cells = lower_guard_cells + 1,
        T = Float64,
    ) where {D,G<:AbstractGrid{D}}

        index_offset = ntuple(_ -> lower_guard_cells + 1, D)
        values = zeros(
            T,
            (grid.num_cells .+ lower_guard_cells .+ upper_guard_cells)...,
            num_elements,
        )

        new{T,num_elements,D,D + 1,G}(
            values,
            grid,
            offset,
            lower_guard_cells,
            upper_guard_cells,
            index_offset,
        )
    end
end

num_elements(f::Field{T,N}) where {T,N} = N

function Base.show(io::IO, f::Field{T,N,D,D1,G}) where {T,N,D,D1,G}
    print(
        io,
        "Field(grid=$(f.grid), offset=$(f.offset),
num_elements=$(N), T=$T",
    )
end

@inline function Base.eachindex(f::Field)
    return CartesianIndices(
        Tuple(UnitRange.(f.index_offset, f.index_offset .+ f.grid.num_cells .- 1)),
    )
end

@inline function Base.getindex(f::Field, I...)
    return f.values[I...]
end

@inline function Base.setindex!(f::Field, v, I...)
    f.values[I...] = v
end

@inline Base.lastindex(f::Field) = lastindex(f.values)
@inline Base.firstindex(f::Field) = firstindex(f.values)

@inline function cell_index_to_cell_coords(f::Field, I)
    # align_offset = f.offset == edge ? 0.5 : 0.0
    # return Tuple(I) .- f.index_offset .+ align_offset
    return Tuple(I) .- f.index_offset
end

@inline function cell_index_to_cell_coords(f::Field{T,N}, I, dim) where {T,N}
    node_cell_coords = cell_index_to_cell_coords(f, I)

    if f.offset == node
        return node_cell_coords
    elseif f.offset == edge
        return node_cell_coords .+ 0.5 .* unit_vec(dim, Val(N))
    elseif f.offset == face
        return node_cell_coords .+ 0.5 .* orth_vec(dim, Val(N))
    end
end

@inline function cell_coords_to_cell_index(f::Field, idxs)
    return CartesianIndex(Tuple(floor.(Int, idxs) .+ f.index_offset))
end

@inline function phys_coords_to_cell_index_ittr(f::Field, xs, width::Int)
    cell_coords = phys_coords_to_cell_coords(f.grid, xs)
    cell_index = Tuple(cell_coords_to_cell_index(f, cell_coords))

    extra_width = 0
    if f.offset == edge || f.offset == face
        extra_width = 1
    end

    # Need to add one to lower bound to account for the fact that the particle
    # is 'in' the cell at cell_index.
    low_bounds = cell_index .- width .+ 1 .- extra_width
    high_bounds = cell_index .+ width .+ extra_width
    return cell_coords, CartesianIndices(Tuple(UnitRange.(low_bounds, high_bounds)))
end

function Base.axes(f::Field)
    return UnitRange.(f.lower_guard_cells .+ 1, f.lower_guard_cells .+ f.grid.num_cells)
end
