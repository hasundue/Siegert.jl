using Test
using Siegert
using Printf: @sprintf
using Logging: @debug

@testset "Phase (optional, square well l=0)" begin
    a = 1.0
    l = 0
    V0 = -112.5
    V = _ -> V0
    N = 19

    ks, = sps_solve(N, l, a, V; b = 0.0, exact_metric = true)
    S_prod = s_from_eigs(ks, a)

    # TON anchor: high-k phase value reported in the paper
    δ_num = 0.5 * atan(imag(S_prod(20.0)), real(S_prod(20.0)))
    # Phase is periodic; shift by integer multiples of π to compare
    δ_ton = -10.681_639_621_0
    m = round(Int, (δ_num - δ_ton) / π)
    δ_adj = δ_num - m * π
    @test isapprox(δ_adj, δ_ton; atol = 1e-6, rtol = 0)

    # Optional low-k spot checks vs analytic S(k)
    ksamp = [0.2, 0.5, 1.0, 1.5]
    for k in ksamp
        S_ana = s_square_well_l0(k, a, V0)
        S_num = S_prod(k)
        @test isapprox(abs(S_num), 1.0; atol = 5e-3, rtol = 0)
        @test isapprox(S_num, S_ana; atol = 2e-2, rtol = 0)
        @debug @sprintf(
            "k=%4.2f: |S_num|=%.6f, |S_ana|=%.6f, |ΔS|=%.2e",
            k,
            abs(S_num),
            abs(S_ana),
            abs(S_num - S_ana)
        )
    end
end
