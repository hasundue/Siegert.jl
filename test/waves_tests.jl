using Test
using Siegert
using LinearAlgebra: norm

@testset "Waves: S via TON (59) has |S|=1 on real k" begin
    a = 1.0;
    l = 0;
    V0 = -112.5;
    V = _ -> V0
    N = 19
    sps = solve_sps(N, l, a, V; b = 0.0)

    S = scattering_matrix(sps)
    for k in (0.2, 0.5, 1.0, 1.5, 2.0)
        @test isapprox(abs(S(k)), 1.0; atol = 1e-10, rtol = 0)
    end
end

@testset "Waves: interior standing wave on SPS grid" begin
    a = 1.0;
    l = 0;
    V0 = -112.5;
    V = _ -> V0
    N = 19
    sps = solve_sps(N, l, a, V; b = 0.0)

    k = 1.0
    u = standing_wave_on_grid(sps, k)
    # Realness on grid (up to FP)
    @test maximum(abs.(imag.(complex.(u)))) < 1e-10

    # Interior residual sanity for outgoing solve (scale by operator norm)
    Qp = sps.H̃ .- (0.0 + im * k * sps.a) .* sps.L .- (k^2) .* sps.ξ
    cout = outgoing_on_grid(sps, k)
    res = Qp * cout
    using LinearAlgebra: opnorm, Symmetric, eigen
    @test norm(res) / (opnorm(Qp) * max(norm(cout), 1e-12)) < 1e-6

    # Boundary constraint v'c ≈ 1 (from augmented LSQ)
    v = ComplexF64.(Siegert._boundary_vector(sps.L))
    @test abs(v' * cout - 1) < 1e-6
end
