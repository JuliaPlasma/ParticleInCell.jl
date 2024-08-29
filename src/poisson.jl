field_solve_three_pt(k, dx) = k^2 * sinc(k * dx / 2 / pi)^2
field_solve_lagrange(k, dx) = k^2 * sinc(k * dx / 2 / pi)^2 * (2 + cos(k * dx)) / 3
field_solve_five_pt(k, dx) = -1 * k^2 * sinc(k * dx / 2 / pi)^2 * (-7 + cos(k * dx)) / 6

struct PoissonSolveFFT{T,D,G,P,F<:AbstractField} <: AbstractSimulationStep
    rho::F
    phi::F

    Ksq_inv::Array{T,D}
    ft_vector::Array{Complex{T},D}

    fft_plan::P

    function PoissonSolveFFT(
        rho::F,
        phi::F,
        k2s = field_solve_three_pt,
    ) where {T,D,G,F<:AbstractField{T,D,G}}
        # This restriction could possibly be relaxed to just require compatible grids...
        @assert rho.grid === phi.grid
        # Currently only supports periodic boundary conditions...
        @assert all(rho.grid.periodic)
        @assert num_elements(rho) == 1 && num_elements(phi) == 1

        epsilon_0 = 8.8541878128e-12
        # Ksq_inv = reshape(copy(rho.values), size(rho.values)[1:end-1]...)
        Ksq_inv = zeros(T, rho.grid.num_cells...)
        ft_vector = zeros(Complex{T}, rho.grid.num_cells...)

        grid = rho.grid
        # TODO: Make this use the cell_lengths logic in Grid
        sim_lengths = grid.upper_bounds .- grid.lower_bounds
        cell_lengths = sim_lengths ./ grid.num_cells
        for I in eachindex(Ksq_inv)
            It = Tuple(I)

            if any(x -> x == 1, It)
                Ksq_inv[I] = 0
                continue
            end

            # ks = 2π .* (It .- 1) ./ sim_lengths
            ks =
                2π .*
                ifelse.(It .<= size(Ksq_inv) ./ 2, It .- 1, It .- size(Ksq_inv) .- 1) ./
                sim_lengths
            inv_Ksqs = (k2s.(ks, cell_lengths)) .^ (-1) ./ epsilon_0

            Ksq_inv[I] = prod(inv_Ksqs)
        end

        fft_plan = plan_fft!(ft_vector)

        new{T,D,G,typeof(fft_plan),F}(rho, phi, Ksq_inv, ft_vector, fft_plan)
    end
end

function step!(step::PoissonSolveFFT)
    step.ft_vector .= view(step.rho.values, eachindex(step.rho))
    step.fft_plan * step.ft_vector
    step.ft_vector .= step.ft_vector .* step.Ksq_inv
    inv(step.fft_plan) * step.ft_vector
    view(step.phi.values, eachindex(step.phi)) .= real.(step.ft_vector)
end
