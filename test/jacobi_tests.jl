using Test
using Siegert
using QuadGK
using LinearAlgebra: Diagonal, Symmetric, diag, isdiag, norm, eigen

# Tolerances
# - Use looser atol for numerical integrations (QuadGK)
# - Use tighter atol for discrete/quadrature identities
const ATOL_INTEGRAL = 1e-4
const ATOL_DISCRETE = 1e-10

@testset "Jacobi L2 basis" begin
    for (α, β) in ((0.0, 0.0), (0.0, 2.0))
        N = 6
        φ = jacobi_basis_L2(N, α, β)
        @test length(φ) == N
        # Orthonormality by numerical integration (QuadGK)
        for m = 1:N, n = 1:N
            g(x) = φ[m](x) * φ[n](x)
            val, _ = quadgk(g, -1.0, 1.0; atol = ATOL_INTEGRAL)
            @test isapprox(val, m == n ? 1.0 : 0.0; atol = ATOL_INTEGRAL, rtol = 0)
        end
    end
end

@testset "Jacobi DVR basis" begin
    for l in (0, 1)
        for N in (3, 4)
            ψ, z, Λ = jacobi_dvr_basis(N, l)
            @test length(ψ) == N
            # Cardinality on Gauss–Jacobi grid: ψ_i(z_j) = δ_ij / √Λ_j
            for i = 1:N, j = 1:N
                val = ψ[i](z[j]) * sqrt(Λ[j])
                @test isapprox(val, i == j ? 1.0 : 0.0; atol = ATOL_DISCRETE, rtol = 0)
            end
            # Orthonormality by numerical integration (QuadGK)
            for i = 1:N, j = 1:N
                g(x) = ψ[i](x) * ψ[j](x)
                val, _ = quadgk(g, -1.0, 1.0; atol = ATOL_INTEGRAL)
                @test isapprox(val, i == j ? 1.0 : 0.0; atol = ATOL_INTEGRAL, rtol = 0)
            end
            # Orthonormality by Gauss–Jacobi quadrature using Λ
            for i = 1:N, j = 1:N
                val = sum(Λ .* (ψ[i].(z) .* ψ[j].(z)))
                @test isapprox(val, i == j ? 1.0 : 0.0; atol = ATOL_DISCRETE, rtol = 0)
            end
        end
    end
end
