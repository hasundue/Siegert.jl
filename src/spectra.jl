# Waves and scattering utilities built on SPS
# Implements scattering observables and interior waves from SPSs.
#
# References
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   "Siegert pseudostate formulation of scattering theory: One-channel case,"
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077
#   (See Sec. IV for scattering observables and interior waves.)
# - D. A. Batishchev and O. I. Tolstikhin,
#   "Siegert pseudostate formulation of scattering theory: Nonzero angular momenta
#   in the one-channel case," Phys. Rev. A 75, 062704 (2007).
#   doi:10.1103/PhysRevA.75.062704 (See Eq. (85) for the product formula.)

using LinearAlgebra: Symmetric, eigen, dot

"""
    scattering_matrix(sps::AbstractSPS) -> (k::Real) -> ComplexF64

Product formula over all `2N + l` SPS momenta. [BT07 Eq.(85)]
"""
function scattering_matrix(sps::AbstractSPS)
    ks = momenta(sps)
    a = box_size(sps)
    return function (k::Real)
        P = prod(((kn + k) / (kn - k) for kn in ks))
        cis(-2k * a) * P
    end
end

function phase_shift(sps::AbstractSPS)
    S = scattering_matrix(sps)
    return function (k::Real)
        z = S(k)
        return 0.5 * atan(imag(z), real(z))
    end
end
