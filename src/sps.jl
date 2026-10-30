# References
# - O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
#   "Siegert pseudostate formulation of scattering theory: One-channel case,"
#   Phys. Rev. A 58, 2077–2095 (1998). doi:10.1103/PhysRevA.58.2077 [TON98]
# - D. A. Batishchev and O. I. Tolstikhin,
#   "Siegert pseudostate formulation of scattering theory: Nonzero angular momenta
#   in the one-channel case," Phys. Rev. A 75, 062704 (2007).
#   doi:10.1103/PhysRevA.75.062704 [BT07]

"""
Map x ∈ [-1, 1] to r ∈ [0, a].
"""
x_to_r(x, a) = a * (1 + x) / 2

"""
Build boundary matrix L = π_i(1)π_j(1) from DVR basis matrix. [TON98 Eq.(C18)]
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
Centrifugal term l(l+1)/2 in the TON98 radial equation. [TON98 Eq.(C4)]

After dividing the TON98 equation by r², this term gives the usual
centrifugal potential l(l+1)/(2r²).
"""
centrifugal_term(l::Real) = l * (l + 1) / 2

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

    U_diag = centrifugal_term(l) .+ (r .^ 2) .* V.(r)
    U = diagonal_dvr(U_diag)

    K̃ = kinetic_dvr(T, φ)
    H̃ = K̃ + U

    return K̃, ρ, L, H̃, dvr, r, π[:, end]
end

"""
Return the zeros z_{l,p}, p = 1..l, of the reverse Bessel polynomial θ_l.

θ_0(z) = 1, θ_1(z) = z + 1, and
θ_l(z) = (2l - 1) θ_{l-1}(z) + z² θ_{l-2}(z) for l ≥ 2. [BT07]
"""
function reverse_bessel_zeros(l::Integer)
    l < 0 && throw(ArgumentError("l must be ≥ 0"))
    l == 0 && return Float64[]
    l == 1 && return [-1.0]

    # Polynomial coefficients in ascending powers of z.
    θ_prev2 = [1.0]      # θ_0
    θ_prev1 = [1.0, 1.0] # θ_1
    θ = θ_prev1
    for n = 2:l
        θ = zeros(Float64, n + 1)
        θ[1:(n)] .+= (2n - 1) .* θ_prev1
        θ[3:(n+1)] .+= θ_prev2
        θ_prev2, θ_prev1 = θ_prev1, θ
    end

    # Monic coefficients in descending powers of z.
    coeffs = reverse(θ)
    C = zeros(Float64, l, l)
    for i = 1:(l-1)
        C[i+1, i] = 1.0
    end
    for i = 1:l
        C[i, end] = -coeffs[end-i+1]
    end

    return ComplexF64.(eigen(C).values)
end

"""
Linearize the QEP to a standard eigenvalue problem. [TON98 Eq.(20,C24); BT07 Eq.(31)]

The angular momentum enters through the zeros z of the reverse Bessel
polynomial and gives a problem of dimension 2N + l.
"""
function sps_linearize_qep(H̃, ρ, L, a, π, z::AbstractVector)
    N = size(H̃, 1)
    size(H̃, 2) == N || throw(ArgumentError("H̃ must be square"))
    size(ρ) == (N, N) || throw(ArgumentError("ρ must be N×N"))
    size(L) == (N, N) || throw(ArgumentError("L must be N×N"))
    length(π) == N || throw(ArgumentError("π must have length N"))

    l = length(z)
    ρ_inv = inv(ρ)
    A = 2 * ρ_inv * (H̃ + L)
    B = -2 * a * ρ_inv * L

    if l == 0
        O = zeros(N, N)
        return [O I; -A -B]
    end

    M = zeros(ComplexF64, 2N + l, 2N + l)
    for i = 1:N
        M[i, N+i] = 1.0
    end
    M[(N+1):(2N), 1:N] .= -A
    M[(N+1):(2N), (N+1):(2N)] .= -B
    coupling = 2 .* (ρ_inv * π)
    for p = 1:l
        M[(N+1):(2N), 2N+p] .= coupling
        M[2N+p, 1:N] .= -(z[p] / a) .* π
        M[2N+p, 2N+p] = -z[p] / a
    end

    return M
end

abstract type AbstractSPS end

"""
SPS (Siegert Pseudostate) solution data structure.

This is the raw eigenproblem result. Average-atom data such as the DVR basis,
normalized coefficients, and normalization constants live in `NormalizedSPS`.

Fields:
- `momenta`: Complex momenta k (Siegert eigenvalues)
- `coefficients`: Raw expansion coefficients C in DVR basis
- `grid`: Grid points r on [0,a]
- `box_size`: Box size a
- `angular_momentum`: Angular momentum quantum number l
- `hamiltonian`: Provisional Hamiltonian H̃ in DVR
- `metric`: Metric matrix ρ in DVR
- `boundary`: Boundary matrix L
"""
struct SPS <: AbstractSPS
    momenta::Vector{ComplexF64}
    coefficients::Matrix{ComplexF64}
    grid::Vector{Float64}
    box_size::Float64
    angular_momentum::Int
    hamiltonian::Matrix{Float64}
    metric::Matrix{Float64}
    boundary::Matrix{Float64}
end

# Accessors shared by raw and average-atom Siegert states.
momenta(s::SPS) = s.momenta
coefficients(s::SPS) = s.coefficients
grid(s::SPS) = s.grid
box_size(s::SPS) = s.box_size
angular_momentum(s::SPS) = s.angular_momentum
hamiltonian(s::SPS) = s.hamiltonian
metric(s::SPS) = s.metric
boundary(s::SPS) = s.boundary

"""
    SPS(N, l, V, a) -> SPS

Solve the SPS quadratic eigenvalue problem for Siegert states. [TON98, BT07]

# Arguments
- `N::Integer`: Number of DVR basis functions
- `l::Integer`: Non-negative angular momentum quantum number
- `V::Function`: Potential energy function V(r)
- `a::Real`: Box size (cutoff radius)

# Returns
- `SPS`: raw Siegert pseudostate solution containing the `2N + l` eigenvalues
  and eigenvectors. Wrap it in `NormalizedSPS` before using the
  average-atom layer.

# Example
```julia
V = square_well(-112.5, 1.0)
sps = SPS(19, 0, V, 1.0)  # s-wave
sps = SPS(19, 1, V, 1.0)  # p-wave
```
"""
function SPS(N::Integer, l::Integer, V::Function, a::Real)
    sps, _ = _solve_sps(N, l, V, a)
    return sps
end

function _solve_sps(N::Integer, l::Integer, V::Function, a::Real)
    l < 0 && throw(ArgumentError("l must be ≥ 0"))
    K̃, ρ, L, H̃, dvr, r, π = sps_matrices(N, l, a, V)
    z = reverse_bessel_zeros(l)
    AL = sps_linearize_qep(H̃, ρ, L, a, π, z)
    F = eigen(AL)
    k = -im .* F.values
    Y = F.vectors
    c = ComplexF64.(Array(@view Y[1:N, :]))
    sps = SPS(
        ComplexF64.(k),
        c,
        Float64.(r),
        Float64(a),
        Int(l),
        Matrix{Float64}(H̃),
        Matrix{Float64}(ρ),
        Matrix{Float64}(L),
    )
    return sps, dvr
end
