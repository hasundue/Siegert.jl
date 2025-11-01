module Siegert

using Jacobi
using LinearAlgebra: Diagonal, eigen, I
using SpecialFunctions: gamma

export jacobi_normsq,
    jacobi_normalization_L2,
    jacobi_basis_L2,
    jacobi_gauss,
    jacobi_weight,
    jacobi_dvr_basis,
    sps_matrices,
    sps_linearize_qep,
    sps_solve,
    square_well_V,
    s_square_well_l0,
    s_from_eigs,
    SPSData,
    solve_sps,
    scattering_matrix,
    phase_shift,
    outgoing_on_grid,
    standing_wave_on_grid

include("jacobi.jl")
include("sps.jl")
include("phase.jl")
include("waves.jl")

end
