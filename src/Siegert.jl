module Siegert

using Jacobi
using SpecialFunctions: gamma

export jacobi_normsq,
    jacobi_normalization_L2, jacobi_basis_L2, jacobi_gauss, jacobi_weight, jacobi_dvr_basis

include("jacobi.jl")

end
