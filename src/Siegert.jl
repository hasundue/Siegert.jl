module Siegert

# References
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   “Siegert pseudostate formulation of scattering theory: One-channel case,”
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077

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
    SPSData,
    solve_sps,
    scattering_matrix,
    phase_shift,
    outgoing_on_grid,
    standing_wave_on_grid

include("jacobi.jl")
include("sps.jl")
include("waves.jl")

end
