# Phase and scattering utilities (kept separate from core SPS)
#
# Provides reference/analytic formulas and helpers to build S(k)
# for simple model potentials, to be used in examples or optional tests.

# S-wave (l = 0) square-well analytic S(k)
# Potential: V(r) = V0 for 0 ≤ r ≤ a, and 0 for r > a
# Units: ℏ = m = 1. Here k is the external wave number.
# Returns: complex S(k) with |S| = 1 for real k ≥ 0.
function s_square_well_l0(k::Real, a::Real, V0::Real)
    # Spherical Bessel j0, n0 and their derivatives
    j0(x) = x == 0 ? 1.0 : sin(x) / x
    j0p(x) = (x * cos(x) - sin(x)) / (x^2)
    n0(x) = x == 0 ? -Inf : -cos(x) / x
    n0p(x) = (sin(x) / x) + (cos(x) / (x^2))

    qa = sqrt(k^2 - 2V0) * a
    ka = k * a
    num = k * j0p(ka) * j0(qa) - sqrt(k^2 - 2V0) * j0(ka) * j0p(qa)
    den = k * n0p(ka) * j0(qa) - sqrt(k^2 - 2V0) * n0(ka) * j0p(qa)
    δ = atan(num, den) # correct quadrant phase shift
    return cis(2δ)     # exp(2iδ)
end

# Build S(k) from SPS eigenvalues via product formula
# With SPS symmetry k ↔ −conj(k), the poles sit at −k_n*:
#   S(k) = exp(-2 i k a) Π_n (k - k_n) / (k + k_n*)
# ks: vector of complex eigenvalues k_n, a: channel radius
s_from_eigs(ks::AbstractVector{<:Complex}, a::Real) = function (k::Real)
    z = one(complex(k))
    for kn in ks
        z *= (k - kn) / (k + conj(kn))
    end
    return z * cis(-2k*a)
end
