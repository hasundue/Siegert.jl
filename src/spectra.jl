# Waves and scattering utilities built on SPS
# Implements scattering observables and interior waves from SPSs.
#
# Reference
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   "Siegert pseudostate formulation of scattering theory: One-channel case,"
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077
#   (See Sec. IV for scattering observables and interior waves.)

using LinearAlgebra: Symmetric, eigen, dot

"""
    scattering_matrix(sps::SPS) -> (k::Real) -> ComplexF64
"""
function scattering_matrix(sps::SPS)
    ks = sps.momenta
    a = sps.box_size
    return function (k::Real)
        P = prod(((kn + k) / (kn - k) for kn in ks))
        cis(-2k * a) * P
    end
end

function phase_shift(sps::SPS)
    S = scattering_matrix(sps)
    return function (k::Real)
        z = S(k)
        return 0.5 * atan(imag(z), real(z))
    end
end
