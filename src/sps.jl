# SPS (Siegert Pseudostates) building blocks
#
# This file provides the grid and multiplicative-operator matrices needed
# for the Tolstikhin–Ostrovsky–Nakamura single-channel SPS method.
#
# Notes
# - Uses the Jacobi Gaussian DVR on x ∈ (-1, 1) with (α, β) = (0, 2l).
# - Maps to the physical radial interval r ∈ (0, a) via r(x) = a (1 + x) / 2.
# - Multiplicative operators (e.g., r², V(r)) are represented diagonally in DVR.
# - The kinetic matrix K̃ is implemented via an FBR→DVR transform using
#   derivatives of the L2-orthonormal Jacobi basis φ_n and Gauss–Jacobi
#   quadrature on the x-domain, with the mapping r(x).
#
# Symbol mapping (paper → code):
# - ρ (metric multiplying k²) → ξ (to avoid overloading ρ)
# - Q(k) = H̃ − (b + i k a) L − k² ξ with H̃ = K̃ + U
# - M = −ξ, C = −im a L, K = H̃ − b L (for linearization)
#
# Public API (initial):
# - sps_matrices(N, l, a, V; exact_metric=false) -> (K̃, ξ, L, H̃, z, Λ, ψ, r)
#   where V is a callable V(r)::Real.
# - sps_linearize_qep(H̃, ξ, L, a; b=0.0) -> (A, B) such that A y = k B y
# - sps_solve(N, l, a, V; b=0.0, exact_metric=false) -> (k, C, meta)

const _ONE = 1.0

# Map x ∈ [-1, 1] to r ∈ [0, a]
_x_to_r(x, a) = a * (1 + x) / 2

# Build the Jacobi DVR grid and basis, then map to r
function _sps_grid(N::Integer, l::Real, a::Real)
    ψ, z, Λ = jacobi_dvr_basis(N, l)
    r = _x_to_r.(z, a)
    return ψ, z, Λ, r
end

# Build boundary matrix L = ψ_i(1) ψ_j(1)
function _boundary_matrix(ψ::Vector{<:Function})
    N = length(ψ)
    v = similar(zeros(N))
    for i = 1:N
        v[i] = ψ[i](1.0)
    end
    return v * v'
end

# Diagonal DVR matrix for a multiplicative operator f(x) (or f(r)) evaluated at grid points
_diagonal_dvr(fvals::AbstractVector) = Diagonal(fvals)

# Centrifugal potential term on r-grid (no singularity since r>0 for Gauss nodes)
_centrifugal_diag(l::Real, r::AbstractVector) = (l * (l + 1)) ./ (2 .* (r .^ 2))

# --- Kinetic energy via FBR→DVR ---

# Derivative of Jacobi polynomial P_n^(α,β)(x) via Jacobi.djacobi
_jacobi_poly_derivative(x::Real, n::Integer, α::Real, β::Real) =
    n == 0 ? 0.0 : Jacobi.djacobi(x, n, α, β)

# Values of φ_n and φ_n' at Gauss–Jacobi nodes for given (N, α, β)
# Returns z, Λ, Φ, dΦ with Φ[n,k] = φ_n(z_k), dΦ[n,k] = φ_n'(z_k)
function _phi_vals_and_derivs(N::Integer, α::Real, β::Real)
    # Gauss–Jacobi nodes/weights and discrete measure for unweighted L2
    z, λ = jacobi_gauss(N, α, β)
    ω = jacobi_weight.(z, Ref(α), Ref(β))
    Λ = λ ./ ω

    # Precompute sqrt weight and its derivative factor at nodes
    sqrtω = sqrt.(ω)
    # d/dx ln sqrt(ω) = -α/(2(1-x)) + β/(2(1+x))
    dlog_sqrtω = @. (-α) / (2 * (1 - z)) + β / (2 * (1 + z))

    # Allocate Φ and dΦ with indices [n, k]
    Φ = zeros(Float64, N, N)
    dΦ = zeros(Float64, N, N)

    # Precompute norms and their inverse square roots
    hs = [jacobi_normsq(n - 1, α, β) for n = 1:N]
    invsqrt_h = 1 ./ sqrt.(hs)

    for n = 1:N
        # Vectorized values across nodes
        Pn = Jacobi.jacobi.(z, n - 1, α, β)
        dPn = n == 1 ? zeros(length(z)) : Jacobi.djacobi.(z, n - 1, α, β)
        scale = invsqrt_h[n] .* sqrtω
        Φ[n, :] = (scale .* Pn)'
        # φ' = sqrtω/√h * (P' + P * d/dx ln sqrtω)
        dΦ[n, :] = (scale .* (dPn .+ Pn .* dlog_sqrtω))'
    end

    return z, Λ, Φ, dΦ
end

# Build kinetic matrix K̃ in the DVR basis using
#
# Derivation of the 1/a factor (from (C15)):
# In a.u., T = −(1/2) d²/dr², and the bilinear form after one integration by parts is
# (1/2) ∫_0^a (du/dr)(dv/dr) dr (boundary term is handled via L).
# With r = a(1+x)/2, dr = (a/2) dx and d/dr = (2/a) d/dx, so
# (1/2)∫ (du/dr)(dv/dr) dr = (1/a) ∫ (du/dx)(dv/dx) dx.
# Discretize in the φ-basis using Gauss–Jacobi quadrature: K_FBR = (1/a) dΦ Diag(Λ) dΦ',
# and map to DVR with S = Diag(√Λ) Φ': K̃ = S K_FBR S'.
function _kinetic_dvr(N::Integer, l::Real, a::Real)
    α = 0.0
    β = 2l
    z, Λ, Φ, dΦ = _phi_vals_and_derivs(N, α, β)

    K_FBR = (1 / a) * (dΦ * Diagonal(Λ) * dΦ')
    S = Diagonal(sqrt.(Λ)) * Φ'
    K̃ = S * K_FBR * S'
    return K̃
end

# Generic multiplicative operator via FBR→DVR using Gauss–Jacobi quadrature.
# Given g(x), build G̃ = S (Φ Diag(Λ) Diag(g(z)) Φ') S', with
# S = Diag(√Λ) Φ'. For g(x) = r(x)^2, this reproduces the DVR diagonal exactly.
function _multiplicative_dvr_from_fbr(N::Integer, l::Real, a::Real, g::Function)
    α = 0.0
    β = 2l
    z, Λ, Φ, _ = _phi_vals_and_derivs(N, α, β)
    W = Diagonal(Λ)
    Dg = Diagonal(g.(z))
    G_FBR = Φ * W * Dg * Φ'
    S = Diagonal(sqrt.(Λ)) * Φ'
    return S * G_FBR * S'
end

# Exact metric option: use FBR→DVR for g(x) = r(x)^2
_metric_exact_dvr(N::Integer, l::Real, a::Real) =
    _multiplicative_dvr_from_fbr(N, l, a, x -> (_x_to_r(x, a))^2)

# Assemble SPS matrices
function sps_matrices(N::Integer, l::Real, a::Real, V::Function; exact_metric::Bool = false)
    N < 1 && throw(ArgumentError("N must be ≥ 1"))
    a <= 0 && throw(ArgumentError("a must be > 0"))

    ψ, z, Λ, r = _sps_grid(N, l, a)

    # Boundary matrix at x = 1
    L = _boundary_matrix(ψ)

    # Metric-like matrix for r^2 factor
    ξ = if exact_metric
        _metric_exact_dvr(N, l, a)
    else
        _diagonal_dvr(r .^ 2)
    end

    # Potential: centrifugal + external V(r), both multiplicative in DVR
    U_diag = _centrifugal_diag(l, r) .+ V.(r)
    U = _diagonal_dvr(U_diag)

    # Kinetic via FBR→DVR
    K̃ = _kinetic_dvr(N, l, a)

    # Provisional Hamiltonian in DVR
    H̃ = K̃ + U

    return K̃, ξ, L, H̃, z, Λ, ψ, r
end

# Linearize QEP: Q(k) c = 0 with Q(k) = k^2 M + k C + K
# Identify M, C, K from (H̃ − (b + i k a) L − k^2 ξ):
#   M = −ξ, C = −im * a * L, K = H̃ − b L
function sps_linearize_qep(H̃, ξ, L, a; b::Real = 0.0)
    N = size(H̃, 1)
    size(H̃, 2) == N || throw(ArgumentError("H̃ must be square"))
    size(ξ) == (N, N) || throw(ArgumentError("ξ must be N×N"))
    size(L) == (N, N) || throw(ArgumentError("L must be N×N"))

    M = -ComplexF64.(ξ)
    C = -im * a * ComplexF64.(L)
    K = ComplexF64.(H̃ .- b .* L)

    Z = zeros(ComplexF64, N, N)
    Ieye = Matrix{ComplexF64}(I, N, N)

    A = [Z Ieye; -K -C]
    B = [Ieye Z; Z M]
    return A, B
end

# Solve the QEP via linearization
function sps_solve(
    N::Integer,
    l::Real,
    a::Real,
    V::Function;
    b::Real = 0.0,
    exact_metric::Bool = false,
)
    K̃, ξ, L, H̃, z, Λ, ψ, r = sps_matrices(N, l, a, V; exact_metric = exact_metric)
    A, B = sps_linearize_qep(H̃, ξ, L, a; b = b)
    F = eigen(A, B)
    k = F.values
    Y = F.vectors
    # Extract c coefficients from lower block of eigenvectors
    C = @view Y[(N+1):end, :]
    return k, C, (; A, B, H̃, ξ, L, z, Λ, ψ, r)
end

# Convenience: square well V(r) = 0 inside [0,a]
const square_well_V = r -> 0.0
