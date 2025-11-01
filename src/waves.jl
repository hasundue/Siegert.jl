# Waves and scattering utilities built on SPS
# Implements scattering observables and interior waves from SPSs.
#
# Reference
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   “Siegert pseudostate formulation of scattering theory: One-channel case,”
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077
#   (See Sec. IV for scattering observables and interior waves.)

using LinearAlgebra: Symmetric, eigen, dot

# Structured SPS result for downstream scattering/waves
struct SPSData
    ks::Vector{ComplexF64}
    C::Matrix{ComplexF64}            # N×(2N)
    r::Vector{Float64}               # length N
    a::Float64
    l::Int
    H̃::Matrix{Float64}              # N×N dense
    ξ::Matrix{Float64}               # N×N dense
    L::Matrix{Float64}               # N×N dense, rank-1
end

"""
    solve_sps(N, l, a, V; b=0.0) -> SPSData

Runs the SPS pipeline and returns a structured result for downstream
scattering and wave reconstruction [TON, Sec. IV].
"""
function solve_sps(N::Integer, l::Integer, a::Real, V::Function; b::Real = 0.0)
    ks, C, meta = sps_solve(N, l, a, V; b = b)
    H̃ = Matrix{Float64}(meta.H̃)
    ρ = Matrix{Float64}(meta.ρ)
    L = Matrix{Float64}(meta.L)
    return SPSData(
        ComplexF64.(ks),
        ComplexF64.(C),
        Float64.(meta.r),
        float(a),
        Int(l),
        H̃,
        ρ,
        L,
    )
end

"""
    scattering_matrix(sps::SPSData; rtol=1e-10) -> (k::Real) -> ComplexF64

Implements TON Eq. (59) directly:
  S(k) = exp(-2 i k a) * (1 + i a k R(k)) / (1 - i a k R(k))
with
  R(k) = Σ_j γ_j^2 / (k - k_j)
where the sum runs over an independent set of Siegert poles and
γ_j = v^T c_j uses the boundary vector v from L ≈ v v'.
"""
function scattering_matrix(sps::SPSData)
    ks = sps.ks
    a = sps.a
    return function (k::Real)
        P = prod(((kn + k) / (kn - k) for kn in ks))
        cis(-2k * a) * P
    end
end

"""
    phase_shift(sps::SPSData; unwrap=false)

- Scalar call: `phase_shift(sps)(k::Real) -> δ(k)::Float64`
- Vector call: `phase_shift(sps)(ks::AbstractVector{<:Real}) -> δs::Vector{Float64}`
  If `unwrap=true`, applies a simple continuous unwrapping across `ks` (assumed sorted).
"""
function phase_shift(sps::SPSData; unwrap::Bool = false)
    S = scattering_matrix(sps)
    scalar = function (k::Real)
        z = S(k)
        return 0.5 * atan(imag(z), real(z))
    end
    function vec(ks::AbstractVector{<:Real})
        δ = [scalar(k) for k in ks]
        if unwrap && !isempty(δ)
            out = similar(δ)
            out[1] = δ[1]
            for i = 2:length(δ)
                d = δ[i] - out[i-1]
                while d > pi/2
                    d -= pi
                end
                while d < -pi/2
                    d += pi
                end
                out[i] = out[i-1] + d
            end
            return out
        else
            return δ
        end
    end
    return k_or_ks -> k_or_ks isa Real ? scalar(k_or_ks) : vec(k_or_ks)
end

# --- Wave reconstruction on the SPS grid ---

# Extract a boundary vector v from L ≈ v v' (rank-1). Scale is arbitrary.
_boundary_vector(L::AbstractMatrix{<:Real}) =
    let
        F = eigen(Symmetric(L))
        i = argmax(F.values)
        sqrt(F.values[i]) * F.vectors[:, i]
    end

# Build Q±(k) = H̃ - (± i k a) L - k^2 ξ (with b=0 convention used elsewhere).
_Q(sps::SPSData, k::Real, sign_im::Int) =
    sps.H̃ .- (0.0 + im * (sign_im * k * sps.a)) .* sps.L .- (k^2) .* sps.ξ

"""
    outgoing_on_grid(sps::SPSData, k::Real; weight=1e3) -> Vector{ComplexF64}

Returns the outgoing interior DVR vector c(k) by solving a constrained
least-squares system:
  minimize ||Q+(k) c|| subject to v' c = 1

Implemented via tall solve `[Q; α v'] c ≈ [0; α]`, where v is the boundary
vector from L and α is `weight`. For moderate k away from poles, this
produces a stable outgoing solution on the SPS grid.
"""
function outgoing_on_grid(sps::SPSData, k::Real; weight::Real = 1e3)
    N = length(sps.r)
    Qp = _Q(sps, k, +1)
    v = ComplexF64.(_boundary_vector(sps.L))
    α = ComplexF64(weight)
    # Augmented least squares: [Q; α v'] c ≈ [0; α]
    M = [Qp; (α .* permutedims(v))]
    b = vcat(zeros(ComplexF64, N), α)
    return M \ b
end

"""
    standing_wave_on_grid(sps::SPSData, k::Real; weight=1e3) -> Vector{Float64}

Forms the real standing wave on the SPS grid via symmetric combination of
outgoing/incoming interior solutions:
  u_stand ≈ 0.5 * (c_out + c_in)
Returns the real part to enforce a real standing wave representation.
"""
function standing_wave_on_grid(sps::SPSData, k::Real; weight::Real = 1e3)
    cout = outgoing_on_grid(sps, k; weight = weight)
    # Incoming: flip the sign on the boundary term
    N = length(sps.r)
    Qm = _Q(sps, k, -1)
    v = ComplexF64.(_boundary_vector(sps.L))
    α = ComplexF64(weight)
    Mm = [Qm; (α .* permutedims(v))]
    bm = vcat(zeros(ComplexF64, N), α)
    cin = Mm \ bm
    u = 0.5 .* (cout .+ cin)
    return Float64.(real.(u))
end
