using Test
using Siegert

@testset "Phase anchor (TON, one-channel square well l=0)" begin
    a = 1.0
    l = 0
    V0 = -112.5
    V = _ -> V0
    N = 19

    sps = solve_sps(N, l, a, V; b = 0.0, exact_metric = true)
    S = scattering_matrix(sps)

    # TON anchor: high-k phase value reported in the paper
    δ_num = 0.5 * atan(imag(S(20.0)), real(S(20.0)))
    δ_ton = -10.681_639_621_0
    m = round(Int, (δ_num - δ_ton) / π)
    δ_adj = δ_num - m * π
    @test isapprox(δ_adj, δ_ton; atol = 1e-6, rtol = 0)

    # Optional low-k spot checks vs analytic S(k) for the square well (l=0)
    # Keep the analytic within the test to avoid public API surface.
    let
        j0(x) = x == 0 ? 1.0 : sin(x) / x
        j0p(x) = (x * cos(x) - sin(x)) / (x^2)
        n0(x) = x == 0 ? -Inf : -cos(x) / x
        n0p(x) = (sin(x) / x) + (cos(x) / (x^2))
        s_square_well_l0(k::Real, a::Real, V0::Real) = begin
            qa = sqrt(k^2 - 2V0) * a
            ka = k * a
            num = k * j0p(ka) * j0(qa) - sqrt(k^2 - 2V0) * j0(ka) * j0p(qa)
            den = k * n0p(ka) * j0(qa) - sqrt(k^2 - 2V0) * n0(ka) * j0p(qa)
            δ = atan(num, den)
            cis(2δ)
        end
        for k in (0.2, 0.5, 1.0, 1.5)
            S_ana = s_square_well_l0(k, a, V0)
            S_num = S(k)
            @test isapprox(abs(S_num), 1.0; atol = 5e-3, rtol = 0)
            @test isapprox(S_num, S_ana; atol = 2e-2, rtol = 0)
        end
    end
end
