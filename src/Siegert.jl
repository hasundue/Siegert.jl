module Siegert

# References
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   “Siegert pseudostate formulation of scattering theory: One-channel case,”
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077

using Jacobi: jacobi, zgrjp, wgrjp
using LinearAlgebra: Diagonal, eigen, I
using SpecialFunctions: gamma

export DVR
export SPS
export square_well
export scattering_matrix, phase_shift

include("dvr.jl")
include("potentials.jl")
include("sps.jl")
include("spectra.jl")

end
