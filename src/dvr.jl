# Reference
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   "Siegert pseudostate formulation of scattering theory: One-channel case,"
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077 [TON98]

"""
DVR basis data structure containing all basis information.

Fields:
- `grid`: Grid points x on [-1,1]
- `basis`: DVR basis matrix π[i,j] (localized basis functions)
- `fbr_basis`: FBR basis matrix φ[n,i] (L2-orthonormal basis functions)
- `transformation`: Transformation matrix T from FBR to DVR
- `coefficients`: Transformation coefficients κ
- `weights`: Gauss quadrature weights ω
"""
struct DVR
    grid::Vector{Float64}
    basis::Matrix{Float64}
    fbr_basis::Matrix{Float64}
    transformation::Matrix{Float64}
    coefficients::Vector{Float64}
    weights::Vector{Float64}
end

"""
Return the Gauss–Radau-Jacobi nodes x, where the node x = 1 is included,
and weights w, for given N, α, β. Returns `(x, w)`.
"""
function jacobi_gauss(N::Integer, α::Real, β::Real)
    N < 1 && throw(DomainError(N, "N must be ≥ 1"))
    x = zgrjp(N, α, β)
    w = wgrjp(x, α, β)
    return x, w
end

"""
Squared L2 norm h_n of P_n^(α,β) on [-1,1] with weight (1-x)^α(1+x)^β.
"""
function jacobi_normsq(n::Integer, α::Real, β::Real)
    n < 0 && throw(DomainError(n, "n must be ≥ 0"))
    num = 2^(α + β + 1) * gamma(n + α + 1) * gamma(n + β + 1)
    den = (2n + α + β + 1) * gamma(n + α + β + 1) * gamma(n + 1)
    return num / den
end

"""
Jacobi weight function ω(x) = (1-x)^α (1+x)^β on [-1,1].
"""
jacobi_weight(x, α, β) = (1 - x)^α * (1 + x)^β

"""
Return N×N matrix of L2-orthonormal Jacobi basis functions φ[n,i] = φ_n(x_i) for n = 0:(N-1). [TON98 Eq.(C6)]
"""
function jacobi_basis(x::Vector{<:Real}, α::Real, β::Real)
    N = length(x)
    w = jacobi_weight.(x, α, β)
    i = 0:(N-1) # n-1 indices for 0-based based functions, being consistent with TON98 (C6)
    h = jacobi_normsq.(i, α, β)
    P = jacobi.(x', i, α, β)
    return sqrt.(w' ./ h) .* P
end

"""
Compute the transformation coefficients κ_i = √(ω_i / w(x_i))
where ω_i are the Gauss–Jacobi weights and w(x_i) is the Jacobi weight at node i.
"""
transformation_coefficients(ω::Vector{<:Real}, w::Vector{<:Real}) = sqrt.(ω ./ w)

"""
Build transformation matrix T from Eq. C10: T_ni = √(λ_i / w(x_i)) ψ_n(x_i)
where ψ_n is the L2-orthonormal basis and w(x_i) is the Jacobi weight at node i.
"""
function transformation_matrix(φ::Matrix{<:Real}, κ::Vector{<:Real})::Matrix{Float64}
    return φ .* κ'
end

"""
Construct the Jacobi-based Gaussian DVR basis {π_i} associated with the Gauss–Jacobi grid for (N, l).
Returns DVRBasis struct containing all basis information.
"""
function DVR(N::Integer, l::Real)
    N < 1 && throw(DomainError(N, "N must be ≥ 1"))
    α = 0
    β = 2l
    x, ω = jacobi_gauss(N, α, β)
    w = jacobi_weight.(x, α, β)
    φ = jacobi_basis(x, α, β)
    κ = transformation_coefficients(ω, w)
    T = transformation_matrix(φ, κ)
    π = T' * φ
    return DVR(x, π, φ, T, κ, ω)
end
