using Test
using Siegert

@testset "Jacobi DVR basis" begin
    for (N, l) in ((4, 0), (5, 1))
        basis = DVR(N, l)
        x, π, κ = basis.grid, basis.basis, basis.coefficients
        @test length(x) == N
        # Cardinality on Gauss–Jacobi grid: π_i(x_j) = δ_ij / κ_j
        for i = 1:N, j = 1:N
            @test isapprox(π[i, j], i == j ? κ[i]^(-1) : 0; atol = 1e-12, rtol = 0)
        end
        # Orthonormality by Gauss–Jacobi quadrature using κ
        for i = 1:N, j = 1:N
            val = sum(κ .^ 2 .* π[i, :] .* π[j, :])
            @test isapprox(val, i == j ? 1.0 : 0.0; atol = 1e-12, rtol = 0)
        end
    end
end
