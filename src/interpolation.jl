struct BSplineChargeInterpolation{S, F, IF} <: AbstractSimulationStep
    species::S
    rho::F

    bspline_order::Int
    interp_width::Int

    interp_func::IF

    function BSplineChargeInterpolation(species::S, rho::F, bspline_order::Int) where {S, F}
        @assert rho.offset == node

        interp_width = div(bspline_order, 2) + 1

        # TODO: calculate b-spline orders in general...
        @assert bspline_order == 1
        interp_func = function(xs)
            weight = one(eltype(xs))
            for x in xs
                if abs(x) > 1
                    return zero(eltype(xs))
                elseif x < 0
                    weight *= 1 + x
                else
                    weight *= 1 - x
                end
            end
            return weight
        end

        new{S, F, typeof(interp_func)}(species, rho, bspline_order,
            interp_width, interp_func)
    end
end

function step!(step::BSplineChargeInterpolation)
    grid = step.rho.grid
    cell_volume = prod(cell_lengths(grid))

    for i in eachindex(step.species.positions)
        # Find which cell the particle is in, and create a CartesianIndices
        # object that extends +/- interp_width in all directions
        Is = phys_coords_to_cell_index_ittr(step.rho,
            step.species.positions[i], step.interp_width)

        # Iterate over nodes within the stencil, and compute the corresponding
        # charge for each node
        for I in Is
            step.rho.values[I] += step.species.charge *
                step.species.weights[i] / cell_volume *
                step.interp_func(
                    interp_dist(step.rho, I, step.species.positions[i])
                )
        end
    end
end

struct BSplineFieldInterpolation{S, F, IF} <: AbstractSimulationStep
    species::S
    elec::F

    bspline_order::Int
    interp_width::Int

    interp_func::IF

    function BSplineFieldInterpolation(species::S, elec::F, bspline_order::Int) where {S, F}
        # TODO: support interpolations from edge/face fields
        @assert elec.offset == node

        interp_width = div(bspline_order, 2) + 1

        # TODO: calculate b-spline orders in general...
        @assert bspline_order == 1
        interp_func = function(xs)
            weight = one(eltype(xs))
            for x in xs
                if abs(x) > 1
                    return zero(eltype(xs))
                elseif x < 0
                    weight *= 1 + x
                else
                    weight *= 1 - x
                end
            end
            return weight
        end

        new{S, F, typeof(interp_func)}(species, elec, bspline_order,
            interp_width, interp_func)
    end
end

function step!(step::BSplineFieldInterpolation)
    grid = step.elec.grid
    cell_sizes = cell_lengths(grid)

    for i in eachindex(step.species.positions)
        # Find which cell the particle is in, and create a CartesianIndices
        # object that extends +/- interp_width in all directions
        Is = phys_coords_to_cell_index_ittr(step.elec,
            step.species.positions[i], step.interp_width)

        # Iterate over nodes within the stencil, and compute the corresponding
        # charge for each node
        for I in Is
            step.species.forces[i] = step.species.forces[i] .+
            step.species.charge .* step.species.weights[i] .*
                step.elec.values[I] .* step.interp_func(
                    interp_dist(step.elec, I, step.species.positions[i])
                )
        end
    end
end
