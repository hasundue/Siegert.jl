# Finite-temperature exchange-correlation models.
#
# References
# - V. V. Karasiev, J. W. Dufty, and S. B. Trickey, "Nonempirical semilocal
#   free-energy density functional for matter under extreme conditions,"
#   Phys. Rev. Lett. 120, 076401 (2018). doi:10.1103/PhysRevLett.120.076401
#   [corrKSDT]. The corrected coefficients are taken from libxc
#   `LDA_XC_CORRKSDT` (ID 318), because the 2018 Supplemental Material was not
#   available; see also V. V. Karasiev, S. B. Trickey, and J. W. Dufty,
#   Phys. Rev. B 99, 195134 (2019) [KSDT-status].
# - Libxc, https://gitlab.com/libxc/libxc, `src/lda_xc_ksdt.c` and
#   `python/lda_exc/lda_xc_ksdt.py`, MPL-2.0.
#
# The unpolarized corrKSDT free energy per particle is
#
#   f_xc(rs, t) = -(aa(t) + bb(t) sqrt(rs) + cc(t) rs) /
#                 (rs (1 + dd(t) sqrt(rs) + ee(t) rs)),
#
# with t = T/T_F, T_F = (3π² n)^(2/3)/2, and rs = (3/(4π n))^(1/3).

const KB_HARTREE_PER_K = 3.166811563455608e-6
const HARTREE_PER_BOHR3_TO_GPA = 2.9421015696522098e4

"""
Abstract supertype for exchange-correlation models.

Concrete models must implement

    xc_potential(xc, n, r, T) -> V_xc

with density `n` and radii `r` on the grid and temperature `T` in Hartree.
They may implement

    xc_free_energy(xc, n, T) -> f_xc

returning the free energy per particle in Hartree. The default returns
`nothing`, which marks a potential-only model and disables
[`thermodynamics`](@ref).
"""
abstract type AbstractXC end

"""
Zero exchange-correlation model. `V_xc = 0` and `f_xc = 0`.
"""
struct NoXC <: AbstractXC end

"""
Corrected KSDT (corrKSDT) unpolarized finite-temperature LDA. See the file
header and [`xc_potential`](@ref) / [`xc_free_energy`](@ref).
"""
struct KSDT <: AbstractXC end

"""
Adapter for a raw callable `xc(n, r) -> V_xc`. The wrapper is potential-only:
`xc_free_energy` returns `nothing`, so it cannot be used with
[`thermodynamics`](@ref).
"""
struct PotentialOnlyXC{F} <: AbstractXC
    potential::F
end

# corrKSDT coefficients from libxc `LDA_XC_CORRKSDT` (unpolarized). The
# corrected 2018 set is used, not the uncorrected KSDT14 set.
const _KSDT_B = (0.342554, 9.141315, 0.448483, 18.553096, 1.05414999729322402165649834296)
const _KSDT_C = (0.875130, -0.256320, 0.953988)
const _KSDT_D = (0.725917, 2.237347, 0.280748, 4.185911, 0.692183)
const _KSDT_E = (0.255415, 0.931933, 0.115398, 17.234117, 0.451437)
const _KSDT_A = (0.750, 3.043630, -0.0922700, 1.703500, 8.310510, 5.11050)
const _KSDT_A0 = 1 / (π * (4 / (9π))^(1 / 3))

# Minimum t used internally. At this point tanh(1/t) is exactly 1 in Float64
# and the exp(-c3/t) term is zero, so this is the t -> 0 limit.
const _KSDT_T_MIN = 1e-300

# Stable sech²(x) = 1 - tanh(x)²; exp(-2|x|) avoids overflow and the
# cancellation in 1 - tanh(x)^2.
@inline function _sech2(x::Real)
    e = exp(-2 * abs(x))
    return 4e / (1 + e)^2
end

# Rational tanh fit used by aa and ee:
#   y = tanh(c/t) * (a1 + a2 t² + a3 t³ + a4 t⁴) /
#                   (1 + a5 t² + a6 t⁴).
@inline function _rational_tanh(t, y, dy)
    p = _KSDT_A[1] + _KSDT_A[2] * t^2 + _KSDT_A[3] * t^3 + _KSDT_A[4] * t^4
    q = 1 + _KSDT_A[5] * t^2 + _KSDT_A[6] * t^4
    dp = 2 * _KSDT_A[2] * t + 3 * _KSDT_A[3] * t^2 + 4 * _KSDT_A[4] * t^3
    dq = 2 * _KSDT_A[5] * t + 4 * _KSDT_A[6] * t^3
    value = y * p / q
    derivative = dy * p / q + y * (dp * q - p * dq) / q^2
    return value, derivative
end

@inline function _rational_tanh(t, y, dy, coeffs)
    p = coeffs[1] + coeffs[2] * t^2 + coeffs[3] * t^4
    q = 1 + coeffs[4] * t^2 + coeffs[5] * t^4
    dp = 2 * coeffs[2] * t + 4 * coeffs[3] * t^3
    dq = 2 * coeffs[4] * t + 4 * coeffs[5] * t^3
    value = y * p / q
    derivative = dy * p / q + y * (dp * q - p * dq) / q^2
    return value, derivative
end

# aa(t) and ee(t): tanh(1/t) times the cubic-in-t² rational function.
@inline function _aa_ee(t, coeffs)
    t = max(t, _KSDT_T_MIN)
    y = tanh(1 / t)
    dy = -_sech2(1 / t) / t^2
    return _rational_tanh(t, y, dy, coeffs)
end

# bb(t) and dd(t): tanh(1/sqrt(t)) times the quartic-in-t² rational function.
@inline function _bb_dd(t, coeffs)
    t = max(t, _KSDT_T_MIN)
    sqrt_t = sqrt(t)
    y = tanh(1 / sqrt_t)
    dy = -_sech2(1 / sqrt_t) / (2 * t * sqrt_t)
    return _rational_tanh(t, y, dy, coeffs)
end

@inline function _aa(t)
    t = max(t, _KSDT_T_MIN)
    y = tanh(1 / t)
    dy = -_sech2(1 / t) / t^2
    value, derivative = _rational_tanh(t, y, dy)
    return _KSDT_A0 * value, _KSDT_A0 * derivative
end

@inline function _ee(t)
    t = max(t, _KSDT_T_MIN)
    y = tanh(1 / t)
    dy = -_sech2(1 / t) / t^2
    return _rational_tanh(t, y, dy, _KSDT_E)
end

@inline function _bb(t)
    return _bb_dd(t, _KSDT_B)
end

@inline function _dd(t)
    return _bb_dd(t, _KSDT_D)
end

@inline function _cc(t)
    t = max(t, _KSDT_T_MIN)
    z = _KSDT_C[3] / t
    e = exp(-z)
    prefactor = _KSDT_C[1] + _KSDT_C[2] * e
    dprefactor = _KSDT_C[2] * z * e / t
    ee, dee = _ee(t)
    return prefactor * ee, dprefactor * ee + prefactor * dee
end

"""
    _corrksdt_f_rs_ft(rs, t) -> (f, f_rs, f_t)

Unpolarized corrKSDT free energy per particle and its partial derivatives at
fixed `t` and fixed `rs`, respectively. `t -> 0` is evaluated through its
stable limit.
"""
@inline function _corrksdt_f_rs_ft(rs::Real, t::Real)
    rs > 0 || throw(DomainError(rs, "rs must be positive"))
    t >= 0 || throw(DomainError(t, "t must be non-negative"))

    aa, daa = _aa(t)
    bb, dbb = _bb(t)
    cc, dcc = _cc(t)
    dd, ddd = _dd(t)
    ee, dee = _ee(t)

    s = sqrt(rs)
    numerator = aa + bb * s + cc * rs
    denominator = rs * (1 + dd * s + ee * rs)
    numerator_rs = bb / (2s) + cc
    numerator_t = daa + dbb * s + dcc * rs
    denominator_rs = 1 + 1.5 * dd * s + 2 * ee * rs
    denominator_t = rs * (ddd * s + dee * rs)

    f = -numerator / denominator
    f_rs = -(numerator_rs * denominator - numerator * denominator_rs) / denominator^2
    f_t = -(numerator_t * denominator - numerator * denominator_t) / denominator^2
    return f, f_rs, f_t
end

"""
    xc_potential(xc::AbstractXC, n, r, T) -> V_xc

Exchange-correlation potential in Hartree on the grid, for density `n`
(bohr^-3), radii `r` (bohr), and temperature `T` (Hartree). The driver applies
the `V_xc[n(R)]` boundary subtraction separately.
"""
function xc_potential(xc::AbstractXC, n, r, T)
    throw(MethodError(xc_potential, (xc, n, r, T)))
end

xc_potential(xc, n, r, T) = xc_potential(_as_xc(xc), n, r, T)

xc_potential(::NoXC, n, r, T) = zeros(Float64, length(n))

function xc_potential(xc::PotentialOnlyXC, n, r, T)
    values = xc.potential(n, r)
    length(values) == length(n) ||
        throw(DimensionMismatch("xc must return one value per grid point"))
    return Float64.(values)
end

function xc_potential(::KSDT, n, r, T)
    length(n) == length(r) || throw(DimensionMismatch("n and r must have the same length"))
    T >= 0 || throw(ArgumentError("T must be non-negative"))
    all(>(0), n) || throw(DomainError(n, "density must be positive"))
    Tout = promote_type(eltype(n), typeof(T), Float64)
    V = Vector{Tout}(undef, length(n))
    for i in eachindex(n)
        rs = (3 / (4π * n[i]))^(1 / 3)
        t_f = (3 * π^2 * n[i])^(2 / 3) / 2
        t = T / t_f
        f, f_rs, f_t = _corrksdt_f_rs_ft(rs, t)
        V[i] = f - (rs / 3) * f_rs - (2 / 3) * t * f_t
    end
    return V
end

"""
    xc_free_energy(xc::AbstractXC, n, T) -> f_xc or nothing

Free energy per particle in Hartree on the grid, or `nothing` for a
potential-only model. Models without a free energy cannot be used with
[`thermodynamics`](@ref).
"""
xc_free_energy(::AbstractXC, n, T) = nothing

xc_free_energy(xc, n, T) = xc_free_energy(_as_xc(xc), n, T)

xc_free_energy(::NoXC, n, T) = zeros(Float64, length(n))

function xc_free_energy(::KSDT, n, T)
    T >= 0 || throw(ArgumentError("T must be non-negative"))
    all(>(0), n) || throw(DomainError(n, "density must be positive"))
    Tout = promote_type(eltype(n), typeof(T), Float64)
    f = Vector{Tout}(undef, length(n))
    for i in eachindex(n)
        rs = (3 / (4π * n[i]))^(1 / 3)
        t_f = (3 * π^2 * n[i])^(2 / 3) / 2
        t = T / t_f
        f[i], _, _ = _corrksdt_f_rs_ft(rs, t)
    end
    return f
end

# Derivatives needed by the thermodynamic observables. Returns
# `(f_xc, f_t, T_F)` or `nothing` when the model has no free energy.
_xc_thermo(::AbstractXC, n, T) = nothing
_xc_thermo(::NoXC, n, T) = (zeros(Float64, length(n)), zeros(Float64, length(n)), nothing)

function _xc_thermo(::KSDT, n, T)
    T >= 0 || throw(ArgumentError("T must be non-negative"))
    all(>(0), n) || throw(DomainError(n, "density must be positive"))
    f = Vector{Float64}(undef, length(n))
    f_t = Vector{Float64}(undef, length(n))
    t_f = Vector{Float64}(undef, length(n))
    for i in eachindex(n)
        rs = (3 / (4π * n[i]))^(1 / 3)
        t_f[i] = (3 * π^2 * n[i])^(2 / 3) / 2
        t = T / t_f[i]
        f[i], _, f_t[i] = _corrksdt_f_rs_ft(rs, t)
    end
    return f, f_t, t_f
end

function _as_xc(xc)
    xc === nothing && return NoXC()
    xc isa AbstractXC && return xc
    xc === :ksdt && return KSDT()
    xc isa Symbol && throw(
        ArgumentError("unknown xc model $xc; use :ksdt, an AbstractXC, or a callable"),
    )
    return PotentialOnlyXC(xc)
end
