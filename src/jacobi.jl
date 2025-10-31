
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

"""
Return the Gauss–Jacobi nodes and weights for given N, α, β.
The weights integrate functions under the measure (1-x)^α (1+x)^β dx.
"""
function jacobi_gauss(N::Integer, α::Real, β::Real)
    N < 1 && throw(DomainError(N, "N must be ≥ 1"))
    z = Jacobi.zgj(N, float(α), float(β))
    w = Jacobi.wgj(z, float(α), float(β))
    return z, w
end

"""
Jacobi weight function ω(x) = (1-x)^α (1+x)^β on [-1,1].
"""
jacobi_weight(x, α, β) = (1 - x)^α * (1 + x)^β

"""
Construct the Jacobi-based Gaussian DVR basis {ψ_i}
associated with the Gauss–Jacobi grid for (N, α, β).

- Each ψ_i(x) lies in span{φ_n}_{n=0}^{N-1} with φ_n orthonormal in L2[-1,1].
- Discrete measure for φ on Gauss–Jacobi grid is Λ_i = λ_i / ω(z_i), where (z, λ) are GJ nodes/weights and ω is the Jacobi weight function.
- Cardinality: ψ_i(z_j) = δ_{ij} / √Λ_j.
- Orthonormality: ∫_{-1}^1 ψ_i(x) ψ_j(x) dx = δ_{ij}.

Returns a vector of length N with callable functions ψ_i.
"""
function jacobi_dvr_basis(N::Integer, α::Real, β::Real)
    N < 1 && throw(DomainError(N, "N must be ≥ 1"))
    φ = jacobi_basis_L2(N, α, β)
    z, λ = jacobi_gauss(N, α, β)
    ω = jacobi_weight.(z, Ref(α), Ref(β))
    Λ = λ ./ ω
    coeffs = [sqrt(Λ[i]) * φ[n](z[i]) for i = 1:N, n = 1:N]
    return [
        let row = view(coeffs, i, :)
            x -> begin
                s = 0.0
                @inbounds for n = 1:N
                    s += row[n] * φ[n](x)
                end
                s
            end
        end for i = 1:N
    ]
end
