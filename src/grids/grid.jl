"""
Parent type for all grid objects, which are used to define the simulation
domain, and to convert between coordinate systems. There are three different
numbering systems that can refer to a location in the simulation domain:
 1. The 'physical coordinates' of a point are the real (dimensionalful)
    coordinates associated with that point. This value can range from the
    lower bounds to the upper bounds of the simulation. This value will
    typically take the form `Vector{T}` or `NTuple{N, T}` where `T <: Real`.
 2. The 'cell coordinates' of a point is the non-dimensional location of the point
    in units of cell lengths. This value can range from 0 to num_cells - 1,
    or outside this range if guard cells are included. The value will
    typically have the type `NTuple{N, Int}`.
 3. The 'cell index' of a point is the `CartesianIndex` that can be used to
    index into field arrays at that point. This value must strictly be
    confined to `axes(field.values)`, which, for any given dimension, will
    typically range from 1 to num_cells + 2*num_guard_cells + 1.

The first two types of indexing, phys_coords and cell_coords, are independent
of the number of guard cells in a given field, and depend only on grid
quantities. Thus utilities for converting between these systems require only a
reference to a grid object. On the other hand, the utilities for cell_index are
specific the field being used, and thus those must be provided an `AbstractField`
to do the coordinate conversion.

In general, physical coordinates are useful when considering the location of a
particle, while the cell index is used to interpolate to and from the particle
locations. The cell coordinates are useful for some interpolation algorithms,
especially those that are defined for non-uniform grids.
"""
abstract type AbstractGrid{D} end

"""
    sim_lengths(grid)

Returns a tuple of the length of the simulation domain in each dimension.
"""
function sim_lengths end

"""
    cell_lengths(grid, [cell_coords])

Returns the length of the cell located at `cell_coords`. For uniform grids, the
`cell_coords` argument is optional.
"""
function cell_lengths end

# Internal helper functions, not to be exported at this time
# Modified from WaterLily.jl
unit_vec(i, ::Val{D}) where {D} = ntuple(j -> j == i ? 1 : 0, D)
orth_vec(i, ::Val{3}) = ntuple(j -> j == i ? 0 : 1, 3)
orth_vec(::Any, ::Union{Val{1},Val{2}}) = error("orth_vec only makes sense in 3D")

"""
    cell_coords_to_phys_coords(grid, idxs, [offset, component])

Converts the cell coordinates `idxs` to a physical coordinate using the geometry
specified in `grid`. Optionally, an offset and component can be specified to get
the physical coordinates of a specific `edge` or `face` of the cell. The component
argument is one-indexed.

For more information on the different types of coordinate systems, see
`AbstractGrid`.
"""
function cell_coords_to_phys_coords end

"""
    phys_coords_to_cell_coords(grid, xs)

Converts the physical coordinate `xs` to a grid coordinate using the geometry
specified in `grid`.

For more information on the different types of coordinate systems, see
`AbstractGrid`.
"""
function phys_coords_to_cell_coords end

include("uniform_cartesian_grid.jl")
