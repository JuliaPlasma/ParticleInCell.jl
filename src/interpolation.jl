struct BSplineChargeInterpolation{S,F,IF} <: AbstractSimulationStep
    species::S
    rho::F

    bspline_order::Int
    interp_width::Int

    interp_func::IF

    function BSplineChargeInterpolation(species::S, rho::F, bspline_order::Int) where {S,F}
        @assert rho.offset == node

        interp_width = div(bspline_order, 2) + 1

        # TODO: calculate b-spline orders in general...
        @assert bspline_order == 1
        interp_func = function (xs)
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

        new{S,F,typeof(interp_func)}(species, rho, bspline_order, interp_width, interp_func)
    end
end

function step!(step::BSplineChargeInterpolation)
    grid = step.rho.grid
    cell_volume = prod(cell_lengths(grid))

    for i in eachindex(step.species.positions)
        # Find which cell the particle is in, and create a CartesianIndices
        # object that extends +/- interp_width in all directions
        Is = phys_coords_to_cell_index_ittr(
            step.rho,
            step.species.positions[i],
            step.interp_width,
        )

        # Iterate over nodes within the stencil, and compute the corresponding
        # charge for each node
        for I in Is
            step.rho.values[I] +=
                step.species.charge * step.species.weights[i] / cell_volume *
                step.interp_func(interp_dist(step.rho, I, step.species.positions[i]))
        end
    end
end

struct BSplineFieldInterpolation{S,F,IF} <: AbstractSimulationStep
    species::S
    elec::F

    bspline_order::Int
    interp_width::Int

    interp_func::IF

    function BSplineFieldInterpolation(species::S, elec::F, bspline_order::Int) where {S,F}
        # TODO: support interpolations from edge/face fields
        @assert elec.offset == node

        interp_width = div(bspline_order, 2) + 1

        # TODO: calculate b-spline orders in general...
        @assert bspline_order == 1
        interp_func = function (xs)
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

        new{S,F,typeof(interp_func)}(
            species,
            elec,
            bspline_order,
            interp_width,
            interp_func,
        )
    end
end

function step!(step::BSplineFieldInterpolation)
    for n in eachindex(step.species.positions)
        # Find which cell the particle is in, and create a CartesianIndices
        # object that extends +/- interp_width in all directions
        Is = phys_coords_to_cell_index_ittr(
            step.elec,
            step.species.positions[n],
            step.interp_width,
        )

        for I in Is
            step.species.forces[n] =
                step.species.forces[n] .+
                step.species.charge .* step.species.weights[n] .* step.elec.values[I] .*
                step.interp_func(interp_dist(step.elec, I, step.species.positions[n]))
        end
    end
end

"""
    compute_knots(degree)

Returns a vector of the knot locations for a polynomial b-spline with `degree`.
"""
function compute_knots(degree)
    collect(range(0, degree + 1)) .- ((degree + 1) / 2)
end

"""
    compute_bspline_coeffs(degree, [T])

Returns a vector of vectors of the coefficients for the b-spline with
polynomial `degree`, and unit-spaced knots, centered at zero.
"""
function compute_bspline_coeffs(degree, T = Float64)
    degree == 0 && return [[one(T)]]

    knots = compute_knots(degree)
    segments = fill(zeros(T, degree + 1), degree + 1)
    lower_segs = compute_bspline_coeffs(degree - 1, T)

    # Special case the two end segments, which only depend on one lower order segment
    segments[1] = poly_integrate(binomial_expand(lower_segs[1], 0.5))
    segments[end] = -1 .* poly_integrate(binomial_expand(lower_segs[end], -0.5))

    # Compute the middle segments as a function of the two neighboring lower order segments
    for i = 2:length(segments)-1
        shifted_splines =
            binomial_expand(lower_segs[i], +0.5) - binomial_expand(lower_segs[i-1], -0.5)
        segments[i] = poly_integrate(shifted_splines)
    end

    # Adjust the constant terms to ensure continuity
    segments[1][1] = 0 - poly_eval(segments[1], knots[1])
    for i = 2:length(segments)
        segments[i][1] =
            poly_eval(segments[i-1], knots[i]) - poly_eval(segments[i], knots[i])
    end

    return segments
end

# Integrate the polynomial specified by coeffs using the standard power
# integration rule. Note that the constant term is undetermined. Thus routine
# sets the constant equal to zero.
function poly_integrate(coeffs)
    new_coeffs = zeros(eltype(coeffs), length(coeffs) + 1)
    for i = 1:length(coeffs)
        new_coeffs[i+1] = coeffs[i] / i
    end
    return new_coeffs
end

function poly_eval(coeffs, x)
    result = zero(x)
    xp = one(x)
    for c in coeffs
        result += c * xp
        xp *= x
    end
    return result
end

# If coeffs = [c_0, c_1, ..., c_n], then create a polynomial that is the
# expansion of sum_n c_n (x + const_term)^(n-1)
function binomial_expand(coeffs, const_term)
    new_coeffs = zeros(eltype(coeffs), length(coeffs))
    for i = 1:length(coeffs), j = 1:i
        new_coeffs[j] += coeffs[i] * binomial(i - 1, j - 1) * const_term^(i - j)
    end
    return new_coeffs
end
