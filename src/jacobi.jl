using Jacobi
using SpecialFunctions: gamma

N = 8 # Demension of the basis
l = 0 # Angular momentum quantum number

# Parameters for Jacobi polynomials
α = 0
β = 2l

"""
Squared L2 norm h_n of P_n^(α,β) on [-1,1] with weight (1-x)^α(1+x)^β.
"""
function jacobi_normsq(n::Integer, α::Real, β::Real)
    n < 0 && throw(DomainError(n, "n must be ≥ 0"))
    a = float(α)
    b = float(β)
    num = 2.0^(a + b + 1) * gamma(n + a + 1) * gamma(n + b + 1)
    den = (2n + a + b + 1) * gamma(n + 1) * gamma(n + a + b + 1)
    return num / den
end

"""
Unweighted normalization factor c_n(x) = √w(x)/√h_n for P_n^(α,β),
so that φ_n(x) = c_n(x) * P_n^(α,β)(x) is orthonormal in L2[-1,1].
"""
function jacobi_normalization_L2(n::Integer, α::Real, β::Real)
    h = jacobi_normsq(n, α, β)
    return x -> sqrt((1 - x)^α * (1 + x)^β) / sqrt(h)
end

"""
Return N L2-orthonormal Jacobi basis functions φ_n(x) for n = 0:(N-1),
with respect to the unweighted inner product on [-1, 1]:
φ_n(x) = √w(x) P_n^(α,β)(x) / √h_n, where w(x)=(1-x)^α(1+x)^β.
"""
function jacobi_basis_L2(N::Integer, α::Real, β::Real)
    N < 0 && throw(DomainError(N, "N must be ≥ 0"))
    return [
        let nn = n
            x -> jacobi_normalization_L2(nn, α, β)(x) * jacobi(x, nn, α, β)
        end for n = 0:(N-1)
    ]
end
