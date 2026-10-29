module Siegert

# References
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   “Siegert pseudostate formulation of scattering theory: One-channel case,”
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077

using Jacobi: zgrjp, wgrjp
using LinearAlgebra: Diagonal, eigen, I
using SpecialFunctions: loggamma

export DVR
export SPS
export square_well
export scattering_matrix, phase_shift
export AADensity
export Contour
export chemical_potential
export common_radial_grid
export contour_parameters
export fermi_dirac
export fixed_potential_density
export matsubara_energies
export normalize_sps
export siegert_density
export siegert_dos
export siegert_green
export siegert_wavefunction
export volume_weights

include("dvr.jl")
include("potentials.jl")
include("sps.jl")
include("spectra.jl")
include("average_atom.jl")

end
