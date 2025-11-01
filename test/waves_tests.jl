using Test
using Siegert
using LinearAlgebra: norm

@testset "Waves: S and phase from SPS (square well)" begin
    a = 1.0;
    l = 0;
    V0 = -112.5;
    V = _ -> V0
    N = 19
    sps = solve_sps(N, l, a, V; b = 0.0, exact_metric = true)

    # S via waves equals product-form
    S1 = scattering_matrix(sps)
    S2 = s_from_eigs(sps.ks, sps.a)
    for k in (0.2, 0.5, 1.0, 1.5)
        @test isapprox(S1(k), S2(k); atol = 5e-12, rtol = 0)
    end

end

@testset "Waves: interior standing wave on SPS grid" begin
    a = 1.0;
    l = 0;
    V0 = -112.5;
    V = _ -> V0
    N = 19
    sps = solve_sps(N, l, a, V; b = 0.0, exact_metric = true)

    k = 1.0
    u = standing_wave_on_grid(sps, k)
    # Realness on grid (up to FP)
    @test maximum(abs.(imag.(complex.(u)))) < 1e-10

    # Interior residual sanity for outgoing solve (scale by operator norm)
    Qp = sps.H̃ .- (0.0 + im * k * sps.a) .* sps.L .- (k^2) .* sps.ξ
    cout = outgoing_on_grid(sps, k)
    res = Qp * cout
    using LinearAlgebra: opnorm, Symmetric, eigen
    @test norm(res) / (opnorm(Qp) * max(norm(cout), 1e-12)) < 1e-3

    # Boundary constraint v'c ≈ 1 (from augmented LSQ)
    v = ComplexF64.(eigen(Symmetric(sps.L)).vectors[:, end])
    @test abs(v' * cout - 1) < 1e-2
end
