using Test
using Siegert
using LinearAlgebra: norm

@testset "Phase shift for square-well potential" begin
    a = 1.0
    l = 0
    V0 = -112.5
    V = square_well(V0, a)
    δ_theory = -10.681_639_621_0
    k = 20.0

    function test_phase_shift(N, rtol)
        sps = SPS(N, l, V, a, -1.0)
        S = scattering_matrix(sps)
        δ_num = 0.5 * atan(imag(S(k)), real(S(k)))
        m = round(Int, (δ_num - δ_theory) / π)
        δ_adj = δ_num - m * π
        @test isapprox(δ_adj, δ_theory; rtol = rtol)
    end

    @testset "2N=38 (6-digit accuracy)" begin
        test_phase_shift(19, 1e-5)
    end

    @testset "2N=50 (12-digit accuracy)" begin
        test_phase_shift(25, 1e-11)
    end
end
