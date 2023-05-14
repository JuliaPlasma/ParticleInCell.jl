struct PoissonSolveFFT{T, D, G, F <: AbstractField} <: AbstractSimulationStep
    rho::F
    phi::F

    Ksq_inv::Array{T, D}
    ft_vector::Array{Complex{T}, D}

    function PoissonSolveFFT(rho::F, phi::F) where {T, D, G, F <: AbstractField{T, D, G}}
        # This restriction could possibly be relaxed to just require compatible grids...
        @assert rho.grid === phi.grid
        # Currently only supports periodic boundary conditions...
        @assert all(rho.grid.periodic)
        @assert num_elements(rho) == 1 && num_elements(phi) == 1

        epsilon_0 = 8.8541878128e-12
        # Ksq_inv = reshape(copy(rho.values), size(rho.values)[1:end-1]...)
        Ksq_inv =   zeros(T, rho.grid.num_cells...)
        ft_vector = zeros(Complex{T}, rho.grid.num_cells...)

        grid = rho.grid
        # TODO: Make this use the cell_lengths logic in Grid
        sim_lengths = grid.upper_bounds .- grid.lower_bounds
        cell_lengths = sim_lengths ./ grid.num_cells
        for I in eachindex(Ksq_inv)
            It = Tuple(I)
            ks = 2π .* (It .- 1) ./ sim_lengths
            grid_angles = ks .* cell_lengths ./ 2
            inv_Ksqs = (cell_lengths ./ (2 .* sin.(grid_angles))).^2 ./ epsilon_0

            Ksq_inv[I] = prod(inv_Ksqs)
        end

        new{T, D, G, F}(rho, phi, Ksq_inv, ft_vector)
    end
end

function step!(::Any, step::PoissonSolveFFT)
    # TODO...
    step.ft_vector .= view(step.rho.values, eachindex(step.rho))
    FFTW.fft!(step.ft_vector)
    step.ft_vector .= step.ft_vector .* step.Ksq_inv
    FFTW.ifft!(step.ft_vector)
    view(step.phi.values, eachindex(step.rho)) .= real.(step.ft_vector)
end

# for i in 1:round(Int, 1 + num_cells / 2)
#     k = 2 * pi * i / sim_length
#     # grid_angle can alternatively be computed as i * pi / num_cells. I
#     # don't do this for intuitions sake.
#     grid_angle = k * cell_length / 2
#     step.ksq_inv[i] = (cell_length / (2 * sin(grid_angle)))^2 / step.epsilon_0
# end

# # Modifiy ft_vector to get phi(k)
# # This mode should be zero if the simulation cell is charge neutral. We
# # assume that it is here.
# step.ft_vector[1] = 0

# # The modes 2:num_cells (of which there are num_cells - 1) are paired in conjugate
# # pairs so rhok[2] = complex_conjugate(rhok[end]) and so on. Thus we need
# # to multiply both by the same scaling factor.
# for i in 1:round(Int, num_cells / 2)
#     k1 = i + 1
#     k2 = num_cells + 1 - i
#     step.ft_vector[k1] = step.ksq_inv[i] * step.ft_vector[k1]
#     k1 == k2 && break
#     step.ft_vector[k2] = step.ksq_inv[i] * step.ft_vector[k2]
# end
