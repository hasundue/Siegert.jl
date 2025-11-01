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
    square_well_V

include("jacobi.jl")
include("sps.jl")

end
