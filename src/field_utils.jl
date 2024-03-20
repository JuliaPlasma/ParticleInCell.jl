struct FiniteDifferenceToEdges{F1,F2,NT} <: AbstractSimulationStep
    nodal_field::F1
    edge_field::F2

    edge_lengths::NT

    function FiniteDifferenceToEdges(
        nodal_field::Field{T1,N1,NodeOffset},
        edge_field::Field{T2,N2,EdgeOffset},
    ) where {T1,T2,N1,N2}
        @assert edge_field.upper_guard_cells >= 1

        edge_lengths = cell_lengths(nodal_field.grid)

        new{typeof(nodal_field),typeof(edge_field),typeof(edge_lengths)}(
            nodal_field,
            edge_field,
            edge_lengths,
        )
    end
end

# TODO: document that this calculates -1*grad(nodal_field). Maybe make it an option?
function step!(step::FiniteDifferenceToEdges{F}) where {T,D,F<:AbstractField{T,D}}
    nodal_field = step.nodal_field
    edge_field = step.edge_field

    # TODO: flip the order of these loops?
    for I in eachindex(edge_field)
        for d = 1:D
            dI = CartesianIndex(unit_vec(d, Val{D}()))
            # TODO: need to update the indexing into edge field, because it has
            # multiple components
            edge_field[I] = -1 * (nodal_field[I+dI] - nodal_field[I]) / step.edge_lengths[d]
        end
    end
end

struct AverageEdgesToNodes{F1,F2} <: AbstractSimulationStep
    edge_field::F1
    nodal_field::F2

    function AverageEdgesToNodes(
        edge_field::Field{T1,N1,EdgeOffset},
        nodal_field::Field{T2,N2,NodeOffset},
    ) where {T1,T2,N1,N2}
        @assert edge_field.lower_guard_cells >= 1
        @assert edge_field.upper_guard_cells >= 1

        new{typeof(edge_field),typeof(nodal_field)}(edge_field, nodal_field)
    end
end

function step!(step::AverageEdgesToNodes{F}) where {T,D,F<:AbstractField{T,D}}
    edge_field = step.edge_field
    nodal_field = step.nodal_field

    for I in eachindex(edge_field)
        for d = 1:D
            dI = CartesianIndex(unit_vec(d, Val{D}()))
            nodal_field[I] = (edge_field[I-dI] + edge_field[I]) / 2
        end
    end
end

struct ZeroField{F} <: AbstractSimulationStep
    field::F
end

step!(step::ZeroField) = step.field.values .= 0
