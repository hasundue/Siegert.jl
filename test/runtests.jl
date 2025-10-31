using Test
using Siegert

# Simple composite Simpson integration on [a,b]
function simpson(f, a, b; n = 2000)
    n % 2 == 0 || (n += 1)
    h = (b - a) / n
    s = f(a) + f(b)
    for i = 1:(n-1)
        x = a + i*h
        s += (i % 2 == 0 ? 2 : 4) * f(x)
    end
    return s * h / 3
end

@testset "Unweighted orthonormal Jacobi basis" begin
    for (α, β) in ((0.0, 0.0), (0.5, 1.0))
        N = 6
        φ = jacobi_basis_L2(N, α, β)
        @test length(φ) == N
        # Check orthonormality
        for m = 1:N, n = 1:N
            g(x) = φ[m](x) * φ[n](x)
            val = simpson(g, -1.0, 1.0)
            @test isapprox(val, m == n ? 1.0 : 0.0; atol = 5e-4, rtol = 5e-4)
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
            # Cardinality: ψ_i(z_j) = δ_ij / √Λ_j
            for i = 1:N, j = 1:N
                val = ψ[i](z[j]) * sqrt(Λ[j])
                @test isapprox(val, i == j ? 1.0 : 0.0; atol = 1e-12, rtol = 1e-12)
            end
            # Orthonormality by Simpson (numerical)
            for i = 1:N, j = 1:N
                g(x) = ψ[i](x) * ψ[j](x)
                val = simpson(g, -1.0, 1.0)
                @test isapprox(val, i == j ? 1.0 : 0.0; atol = 5e-4, rtol = 5e-4)
            end
            # Orthonormality by Gauss–Jacobi quadrature using Λ
            for i = 1:N, j = 1:N
                val = sum(Λ .* (ψ[i].(z) .* ψ[j].(z)))
                @test isapprox(val, i == j ? 1.0 : 0.0; atol = 1e-12, rtol = 1e-12)
            end
        end
    end
end
