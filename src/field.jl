abstract type AbstractField{T,N,O,D,G} end

abstract type GridOffset end
struct NodeOffset <: GridOffset end
struct EdgeOffset <: GridOffset end
struct FaceOffset <: GridOffset end
struct CenterOffset <: GridOffset end

"""
    Field(grid, offset, num_elements,[ lower_guard_cells = 0,
        [upper_guard_cells = lower_guard_cells + 1]])

Stores the values of a field.
"""
struct Field{T,N,O,D,D1,G} <: AbstractField{T,N,O,D,G}
    values::Array{T,D1}
    grid::G
    offset::O

    lower_guard_cells::Int
    upper_guard_cells::Int
    index_offset::NTuple{D,Int}

    function Field(
        grid::G,
        offset::O,
        num_elements::Int,
        lower_guard_cells = 0,
        upper_guard_cells = lower_guard_cells + 1,
        T = Float64,
    ) where {D,O<:GridOffset,G<:AbstractGrid{D}}

        index_offset = ntuple(_ -> lower_guard_cells + 1, D)
        values = zeros(
            T,
            (grid.num_cells .+ lower_guard_cells .+ upper_guard_cells)...,
            num_elements,
        )

        new{T,num_elements,O,D,D + 1,G}(
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

function Base.show(io::IO, f::Field{T,N,O,D,D1,G}) where {T,N,O,D,D1,G}
    print(
        io,
        "Field(grid=$(f.grid), offset=$(O),
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
    return Tuple(I) .- f.index_offset
end

@inline function cell_index_to_cell_coords(f::Field{T,N,NodeOffset}, I, dim) where {T,N}
    return cell_index_to_cell_coords(f, I)
end

@inline function cell_index_to_cell_coords(f::Field{T,N,EdgeOffset}, I, dim) where {T,N}
    return cell_index_to_cell_coords(f, I) .+ 0.5 .* unit_vec(dim, Val(N))
end

@inline function cell_index_to_cell_coords(f::Field{T,N,FaceOffset}, I, dim) where {T,N}
    return cell_index_to_cell_coords(f, I) .+ 0.5 .* orth_vec(dim, Val(N))
end

@inline function cell_index_to_cell_coords(f::Field{T,N,CenterOffset}, I, dim) where {T,N}
    return cell_index_to_cell_coords(f, I) .+ 0.5
end

@inline function cell_coords_to_cell_index(f::Field, idxs)
    return CartesianIndex(Tuple(floor.(Int, idxs) .+ f.index_offset))
end

extra_width_func(f::Field{T,N,NodeOffset}) where {T,N} = 0
extra_width_func(f) = 1

@inline function phys_coords_to_cell_index_ittr(f::Field, xs, width::Int)
    cell_coords = phys_coords_to_cell_coords(f.grid, xs)
    cell_index = Tuple(cell_coords_to_cell_index(f, cell_coords))

    extra_width = extra_width_func(f)

    # Need to add one to lower bound to account for the fact that the particle
    # is 'in' the cell at cell_index.
    low_bounds = cell_index .- width .+ 1 .- extra_width
    high_bounds = cell_index .+ width .+ extra_width
    return cell_coords, CartesianIndices(Tuple(UnitRange.(low_bounds, high_bounds)))
end

function Base.axes(f::Field)
    return UnitRange.(f.lower_guard_cells .+ 1, f.lower_guard_cells .+ f.grid.num_cells)
end
