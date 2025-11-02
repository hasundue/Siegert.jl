# SPS (Siegert Pseudostates) building blocks
#
# This file provides the grid and multiplicative-operator matrices needed
# for the Tolstikhin–Ostrovsky–Nakamura single-channel SPS method.
#
# Reference
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   "Siegert pseudostate formulation of scattering theory: One-channel case,"
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077
#
# Notes
# - Uses the Jacobi Gaussian DVR on x ∈ (-1, 1) with (α, β) = (0, 2l).
# - Maps to the physical radial interval r ∈ (0, a) via r(x) = a (1 + x) / 2.
# - Both kinetic K̃ and metric ρ matrices are implemented via FBR→DVR transform
#   using the TON 1998 Appendix C formalism (Eqs. C10, C20, C21).
# - Multiplicative potential operators V(r) are diagonal in DVR.
#
# Symbol mapping (paper → code):
# - ρ (metric multiplying k²/2) in paper is ρ in code
# - Eq. (C15): QEP (H̃ − (b + ika)L − k²/2 ρ)c = 0 with H̃ = K̃ + U
# - Following TON 1998 Eq. (20) with matrices from Eq. (C24)
#
# Public API (initial):
# - sps_matrices(N, l, a, V) -> (K̃, ρ, L, H̃, z, Λ, ψ, r) where V is V(r)::Real
# - sps_linearize_qep(H̃, ρ, L, a; b=0.0) -> Ã such that Ã y = λ y
# - sps_solve(N, l, a, V; b=0.0) -> (k, C, meta)

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

# Build transformation matrix T from Eq. C10: T_ni = √(λ_i / w(x_i)) ψ_n(x_i)
# where ψ_n is the L2-orthonormal basis and w(x_i) is the Jacobi weight at node i.
function _transformation_matrix_T(N::Integer, α::Real, β::Real)
    z, λ = jacobi_gauss(N, α, β)
    ω = jacobi_weight.(z, Ref(α), Ref(β))
    Λ = λ ./ ω
    φ_vals = jacobi_basis_at_1(N, α, β)  # Will recompute properly below

    # Get φ values at all nodes
    φ = jacobi_basis_L2(N, α, β)
    Φ = zeros(N, N)
    for n = 1:N
        for i = 1:N
            Φ[n, i] = φ[n](z[i])
        end
    end

    # T[n, i] = sqrt(Λ[i]) * φ_n(z[i])
    T = zeros(N, N)
    for n = 1:N
        for i = 1:N
            T[n, i] = sqrt(Λ[i]) * Φ[n, i]
        end
    end

    return T
end

function _kinetic_fbr_ton(N::Integer, α::Real, β::Real, a::Real)
    # Build kinetic matrix K̃^(ψ) in FBR via Eq. (C21a) and (C21b)
    #
    # Reference: Tolstikhin et al., Phys. Rev. A 58, 2077 (1998), Appendix C
    # - Eq. (C21a): K̃^(ψ)_nm = φ_n(1)φ_m(1)[2Σ_{k=1}^{n-1} φ_k²(1) + φ_n²(1) - 1/2] for n < m
    # - Eq. (C21b): K̃^(ψ)_nn = 2φ_n²(1)Σ_{k=1}^{n-1} φ_k²(1) + 1/2(φ_n²(1) - 1/2)²
    # - By symmetry: K̃^(ψ)_nm = K̃^(ψ)_mn for n > m

    # Get boundary values φ_n(1) (L2-orthonormal basis)
    φ_at_1 = jacobi_basis_at_1(N, α, β)

    # Build K̃^(ψ) matrix
    K_psi = zeros(N, N)

    for n = 1:N
        # Compute cumulative sum Σ_{k=1}^{n-1} φ_k²(1)
        cum_sum = n == 1 ? 0.0 : sum(φ_at_1[k]^2 for k = 1:(n-1))

        # Diagonal element: Eq. (C21b)
        K_psi[n, n] = 2 * φ_at_1[n]^2 * cum_sum + 0.5 * (φ_at_1[n]^2 - 0.5)^2

        # Off-diagonal elements: Eq. (C21a) for n < m, then use symmetry
        for m = (n+1):N
            K_psi[n, m] = φ_at_1[n] * φ_at_1[m] * (2 * cum_sum + φ_at_1[n]^2 - 0.5)
            K_psi[m, n] = K_psi[n, m]  # Symmetry
        end
    end

    return K_psi
end

function _kinetic_dvr(N::Integer, l::Real, a::Real)
    α = 0.0
    β = 2l

    # Build transformation matrix T (Eq. C10)
    T = _transformation_matrix_T(N, α, β)

    # Build K̃^(ψ) in FBR (Eq. C21)
    K_psi = _kinetic_fbr_ton(N, α, β, a)

    # Transform to DVR (Eq. C20)
    K̃ = T * K_psi * T

    return K̃
end

# Build metric matrix ρ in the DVR basis using TON 1998 Appendix C method
# ρ_ij = Σ_{n,m} T_ni ρ^(ψ)_nm T_mj  (Eq. C20 structure)
#
# Reference: Tolstikhin et al., Phys. Rev. A 58, 2077 (1998), Appendix C
# - The metric ρ^(ψ) in FBR is computed analogously to K̃^(ψ)
# - For the coordinate r², we have ρ^(ψ)_nm = ∫ ψ_n(x) r(x)² ψ_m(x) dx
# - Using r(x) = a(1+x)/2, this gives ρ^(ψ)_nm = a² ∫ ψ_n(x) (1+x)²/4 ψ_m(x) dx
function _metric_dvr(N::Integer, l::Real, a::Real)
    α = 0.0
    β = 2l

    # Build transformation matrix T (Eq. C10)
    T = _transformation_matrix_T(N, α, β)

    # Get φ values at all nodes for computing FBR matrix elements
    z, λ = jacobi_gauss(N, α, β)
    φ = jacobi_basis_L2(N, α, β)

    # Build ρ^(ψ) in FBR: ρ^(ψ)_nm = a²/4 ∫ ψ_n(x) (1+x)² ψ_m(x) dx
    # Using Gauss quadrature with measure Λ = λ/w (unweighted L2 measure)
    ω = jacobi_weight.(z, Ref(α), Ref(β))
    Λ = λ ./ ω

    # (1+x)² evaluated at quadrature points
    r_sq_vals = [(a * (1 + xi) / 2)^2 for xi in z]

    # Matrix element: ∫ φ_n(x) r(x)² φ_m(x) dx using quadrature
    ρ_psi = zeros(N, N)
    for n = 1:N
        for m = 1:N
            ρ_psi[n, m] = sum(φ[n](z[i]) * r_sq_vals[i] * φ[m](z[i]) * Λ[i] for i = 1:N)
        end
    end

    # Transform to DVR (analogous to Eq. C20)
    ρ = T' * ρ_psi * T

    return ρ
end

# Assemble SPS matrices
function sps_matrices(N::Integer, l::Real, a::Real, V::Function)
    N < 1 && throw(ArgumentError("N must be ≥ 1"))
    a <= 0 && throw(ArgumentError("a must be > 0"))

    ψ, z, Λ, r = _sps_grid(N, l, a)

    # Boundary matrix at x = 1
    L = _boundary_matrix(ψ)

    # Metric via FBR→DVR
    ρ = _metric_dvr(N, l, a)

    # Potential: centrifugal + external V(r), both multiplicative in DVR
    U_diag = _centrifugal_diag(l, r) .+ V.(r)
    U = _diagonal_dvr(U_diag)

    # Kinetic via FBR→DVR
    K̃ = _kinetic_dvr(N, l, a)

    # Provisional Hamiltonian in DVR
    H̃ = K̃ + U

    return K̃, ρ, L, H̃, z, Λ, ψ, r
end

# Linearize QEP following TON 1998 Eq. (20) with matrices from Eq. (C24)
#
# Reference: Tolstikhin et al., Phys. Rev. A 58, 2077 (1998)
# - Eq. (C15): QEP (H̃ − (b + ika)L − k²/2 ρ)c = 0
# - Multiply by 2ρ⁻¹: (2ρ⁻¹(H̃ − bL) − 2ikaρ⁻¹L − k²I)c = 0
# - Standard form: (A + λB + λ²I)c = 0 with λ = ik
# - Eq. (C24): A = 2ρ⁻¹(H̃ − bL), B = −2aρ⁻¹L (factor 2 accounts for 1/2 in C15)
# - Eq. (20): Companion linearization [0 I; -K -C][c; λc] = λ[c; λc]
# - Standard companion form: K = −A, C = −B
# - Solve eigenvalue problem: Ã y = λ y with y = [c; λc]
#
# Returns Ã for the linearized system where eigenvalue λ = ik
function sps_linearize_qep(H̃, ρ, L, a; b::Real = 0.0)
    N = size(H̃, 1)
    size(H̃, 2) == N || throw(ArgumentError("H̃ must be square"))
    size(ρ) == (N, N) || throw(ArgumentError("ρ must be N×N"))
    size(L) == (N, N) || throw(ArgumentError("L must be N×N"))

    # Eq. (C24): Matrices A and B (non-symmetric version)
    ρ_inv = inv(ρ)
    A = 2 * ρ_inv * (H̃ - b * L)
    B = -2 * a * ρ_inv * L

    # Eq. (20): Standard companion linearization
    # For (A + λB + λ²I)c = 0, the companion form is:
    # [0   I][ c]    [ c]
    # [-A -B][λc] = λ[λc]
    Z = zeros(ComplexF64, N, N)
    # Ieye = Matrix{ComplexF64}(I, N, N)

    # K = -ComplexF64.(A)
    # C = -ComplexF64.(B)

    # Ã y = λ y
    Ã = [Z I; -A -B]
    return Ã
end

# Solve the QEP via linearization
function sps_solve(N::Integer, l::Real, a::Real, V::Function; b::Real = 0.0)
    K̃, ρ, L, H̃, z, Λ, ψ, r = sps_matrices(N, l, a, V)
    Ã = sps_linearize_qep(H̃, ρ, L, a; b = b)
    F = eigen(Ã)
    # Swap imag and real to get k from λ = ik
    k = -im .* F.values
    Y = F.vectors
    c = Array(@view Y[1:N, :])
    return k, c, (; Ã, H̃, ρ, L, z, Λ, ψ, r)
end

# Convenience: square well V(r) = 0 inside [0,a]
const square_well_V = r -> 0.0
