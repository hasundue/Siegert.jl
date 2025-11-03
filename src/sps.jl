# Reference
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   "Siegert pseudostate formulation of scattering theory: One-channel case,"
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077 [TON98]

"""
Map x ∈ [-1, 1] to r ∈ [0, a].
"""
x_to_r(x, a) = a * (1 + x) / 2

"""
Build boundary matrix L = π_i(1)π_j(1) from DVR basis matrix. [TON98 Eq.(C19)]
"""
function boundary_matrix(π::Matrix{<:Real})
    π_at_1 = π[:, end]
    return π_at_1 * π_at_1'
end

"""
Diagonal DVR matrix for multiplicative operator f(r) at grid points.
"""
diagonal_dvr(fvals::AbstractVector) = Diagonal(fvals)

"""
Centrifugal potential l(l+1)/(2r²) on r-grid. [TON98 Eq.(1)]
"""
centrifugal_potential(l::Real, r::AbstractVector) = (l * (l + 1)) ./ (2 .* r .^ 2)

"""
Build kinetic matrix K̃^(φ) in FBR using φ basis. [TON98 Eq.(C21a,C21b)]
"""
function kinetic_fbr(φ::Matrix{<:Real})
    φ_at_1 = φ[:, end]
    N = size(φ, 1)
    K_φ = zeros(N, N)

    for n = 1:N
        cum_sum = n == 1 ? 0.0 : sum(φ_at_1[k]^2 for k = 1:(n-1))

        K_φ[n, n] = 2 * φ_at_1[n]^2 * cum_sum + 0.5 * (φ_at_1[n]^2 - 0.5)^2

        for m = (n+1):N
            K_φ[n, m] = φ_at_1[n] * φ_at_1[m] * (2 * cum_sum + φ_at_1[n]^2 - 0.5)
            K_φ[m, n] = K_φ[n, m]
        end
    end

    return K_φ
end

"""
Build kinetic matrix K̃ in DVR via FBR→DVR transformation. [TON98 Eq.(C20)]
"""
function kinetic_dvr(T::Matrix{<:Real}, φ::Matrix{<:Real})
    K_φ = kinetic_fbr(φ)
    K̃ = T' * K_φ * T
    return K̃
end

"""
Build metric matrix ρ in DVR. [TON98 Eq.(C22,C23)]
"""
function metric_dvr(x::Vector{<:Real}, T::Matrix{<:Real}, l::Real, a::Real)
    N = length(x)
    β = 2l

    ρ = Diagonal((1 .+ x) .^ 2)

    Δ_N = 4 * N^2 * (N + β)^2 / ((2N + β)^2 * ((2N + β)^2 - 1))

    T_N = T[N, :]

    ρ += Δ_N * (T_N * T_N')
    ρ *= (a^2) / 4

    return Matrix(ρ)
end

"""
Assemble SPS matrices for TON single-channel method.
"""
function sps_matrices(N::Integer, l::Real, a::Real, V::Function)
    N < 1 && throw(ArgumentError("N must be ≥ 1"))
    a <= 0 && throw(ArgumentError("a must be > 0"))

    dvr = DVR(N, l)
    x, π, φ, T = dvr.grid, dvr.basis, dvr.fbr_basis, dvr.transformation
    r = x_to_r.(x, a)

    L = boundary_matrix(π)
    ρ = metric_dvr(x, T, l, a)

    U_diag = centrifugal_potential(l, r) .+ (r .^ 2) .* V.(r)
    U = diagonal_dvr(U_diag)

    K̃ = kinetic_dvr(T, φ)
    H̃ = K̃ + U

    return K̃, ρ, L, H̃, x, r
end

"""
Linearize QEP to standard eigenvalue problem. [TON98 Eq.(20,C24)]
"""
function sps_linearize_qep(H̃, ρ, L, a; b::Real = 0.0)
    N = size(H̃, 1)
    size(H̃, 2) == N || throw(ArgumentError("H̃ must be square"))
    size(ρ) == (N, N) || throw(ArgumentError("ρ must be N×N"))
    size(L) == (N, N) || throw(ArgumentError("L must be N×N"))

    ρ_inv = inv(ρ)
    A = 2 * ρ_inv * (H̃ - b * L)
    B = -2 * a * ρ_inv * L
    O = zeros(N, N)

    return [O I; -A -B]
end

"""
SPS (Siegert Pseudostate) solution data structure.

Fields:
- `momenta`: Complex momenta k (Siegert eigenvalues)
- `coefficients`: Expansion coefficients C in DVR basis
- `grid`: Grid points r on [0,a]
- `box_size`: Box size a
- `angular_momentum`: Angular momentum quantum number l
- `hamiltonian`: Provisional Hamiltonian H̃ in DVR
- `metric`: Metric matrix ρ in DVR
- `boundary`: Boundary matrix L
"""
struct SPS
    momenta::Vector{ComplexF64}
    coefficients::Matrix{ComplexF64}
    grid::Vector{Float64}
    box_size::Float64
    angular_momentum::Int
    hamiltonian::Matrix{Float64}
    metric::Matrix{Float64}
    boundary::Matrix{Float64}
end

"""
    SPS(N, l, V, a, b) -> SPS

Solve SPS quadratic eigenvalue problem for Siegert states. [TON98]

# Arguments
- `N::Integer`: Number of DVR basis functions
- `l::Integer`: Angular momentum quantum number
- `V::Function`: Potential energy function V(r)
- `a::Real`: Box size (cutoff radius)
- `b::Real`: Boundary parameter (-1 for 3D spherical waves, see TON98 Appendix C, Eq. C1c)

# Returns
- `SPS`: Siegert pseudostate solution containing eigenvalues and eigenvectors

# Example
```julia
V = square_well(-112.5, 1.0)
sps = SPS(19, 0, V, 1.0, -1.0)  # 3D spherical wave with s-wave (l=0)
```
"""
function SPS(N::Integer, l::Integer, V::Function, a::Real, b::Real)
    K̃, ρ, L, H̃, x, r = sps_matrices(N, l, a, V)
    AL = sps_linearize_qep(H̃, ρ, L, a; b = b)
    F = eigen(AL)
    k = -im .* F.values
    Y = F.vectors
    c = Array(@view Y[1:N, :])
    return SPS(
        ComplexF64.(k),
        ComplexF64.(c),
        Float64.(r),
        Float64(a),
        Int(l),
        Matrix{Float64}(H̃),
        Matrix{Float64}(ρ),
        Matrix{Float64}(L),
    )
end
