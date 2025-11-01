using Test
using Siegert
using QuadGK
using LinearAlgebra: Diagonal, Symmetric, diag, isdiag, norm, eigen

atol = 1e-4

@testset "Jacobi L2 basis" begin
    for (α, β) in ((0.0, 0.0), (0.0, 2.0))
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
    for l in (0, 1)
        for N in (3, 4)
            ψ, z, Λ = jacobi_dvr_basis(N, l)
            @test length(ψ) == N
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

@testset "SPS matrices (grid, L, ξ, H̃, K̃)" begin
    for (N, l, a) in ((4, 0, 2.0), (5, 1, 3.0))
        K̃, ξ, L, H̃, z, Λ, ψ, r = sps_matrices(N, l, a, square_well_V)
        # Basic shapes and monotonic grid
        @test size(K̃) == (N, N)
        @test size(ξ) == (N, N) && isdiag(ξ)
        @test size(L) == (N, N)
        @test size(H̃) == (N, N)
        @test length(z) == N && length(Λ) == N && length(ψ) == N && length(r) == N
        @test all(-1 .< z .< 1)
        @test all(Λ .> 0)
        @test all(0 .< r .< a)
        @test issorted(r)

        # Boundary matrix equals ψ(1) ψ(1)'
        v = [ψ[i](1.0) for i = 1:N]
        @test isapprox(L, v * v'; atol = 1e-10, rtol = 0)

        # ξ is r^2 on the diagonal (within FP tolerance)
        @test isapprox(diag(ξ), r .^ 2; atol = 1e-12, rtol = 0)

        # K̃ should be symmetric positive semidefinite
        @test isapprox(K̃, K̃'; atol = 1e-12, rtol = 0)
        vals = eigen(Symmetric(K̃)).values
        @test all(vals .>= -1e-12)
    end

    # Check 1/a scaling of K̃ numerically for fixed (N,l)
    let N = 5, l = 0
        K2, = begin
            K̃, = sps_matrices(N, l, 2.0, square_well_V)
        end
        K3, = begin
            K̃, = sps_matrices(N, l, 3.0, square_well_V)
        end
        # Expect K̃(a) ∝ 1/a: 2*K̃(2) ≈ 3*K̃(3) in Frobenius norm
        A = 2 .* K2
        B = 3 .* K3
        relerr = norm(A .- B) / max(norm(B), 1e-12)
        @test relerr < 5e-2
    end
end

@testset "QEP linearization (residual check)" begin
    N, l, a = 4, 0, 2.0
    K̃, ξ, L, H̃, z, Λ, ψ, r = sps_matrices(N, l, a, square_well_V)
    A, B = sps_linearize_qep(H̃, ξ, L, a; b = 0.0)
    F = eigen(A, B)
    k = F.values
    Y = F.vectors
    # Build residuals Q(k) c with Q(k) = H̃ - (0 + im*k*a)L - k^2 ξ
    res = zeros(length(k))
    for j = 1:length(k)
        kj = k[j]
        c = Y[(N+1):end, j]
        rj = (H̃ - (0.0 + im*kj*a) .* L - (kj^2) .* ξ) * c
        res[j] = norm(rj)
    end
    @test all(res .< 1e-8)
end
