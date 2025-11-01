using Test
using Siegert
using QuadGK

atol = 1e-4

@testset "Unweighted orthonormal Jacobi basis" begin
    for (α, β) in ((0.0, 0.0), (0.5, 1.0))
        N = 6
        φ = jacobi_basis_L2(N, α, β)
        @test length(φ) == N
        # Orthonormality by numerical integration (QuadGK)
        for m = 1:N, n = 1:N
            g(x) = φ[m](x) * φ[n](x)
            val, _ = quadgk(g, -1.0, 1.0; atol = atol)
            @test isapprox(val, m == n ? 1.0 : 0.0; atol = atol, rtol = 0)
        end
    end
end

@testset "Jacobi DVR basis" begin
    for (α, β) in ((0.0, 0.0), (0.5, 1.0))
        for N in (3, 5, 8)
            ψ = jacobi_dvr_basis(N, α, β)
            @test length(ψ) == N
            z, λ = jacobi_gauss(N, α, β)
            ω = jacobi_weight.(z, Ref(α), Ref(β))
            Λ = λ ./ ω
            # Cardinality on Gauss–Jacobi grid: ψ_i(z_j) = δ_ij / √Λ_j
            for i = 1:N, j = 1:N
                val = ψ[i](z[j]) * sqrt(Λ[j])
                @test isapprox(val, i == j ? 1.0 : 0.0; atol = atol, rtol = 0)
            end
            # Orthonormality by numerical integration (QuadGK)
            for i = 1:N, j = 1:N
                g(x) = ψ[i](x) * ψ[j](x)
                val, _ = quadgk(g, -1.0, 1.0; atol = atol)
                @test isapprox(val, i == j ? 1.0 : 0.0; atol = atol, rtol = 0)
            end
            # Orthonormality by Gauss–Jacobi quadrature using Λ
            for i = 1:N, j = 1:N
                val = sum(Λ .* (ψ[i].(z) .* ψ[j].(z)))
                @test isapprox(val, i == j ? 1.0 : 0.0; atol = atol, rtol = 0)
            end
        end
    end
end
