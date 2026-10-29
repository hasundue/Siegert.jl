using Test
using Siegert
using SpecialFunctions: besselj, bessely

spherical_j(l, x) = sqrt(π / (2x)) * besselj(l + 0.5, x)
spherical_y(l, x) = sqrt(π / (2x)) * bessely(l + 0.5, x)

function spherical_j_derivative(l, x)
    l == 0 && return -spherical_j(1, x)
    return spherical_j(l - 1, x) - (l + 1) / x * spherical_j(l, x)
end

function spherical_y_derivative(l, x)
    l == 0 && return -spherical_y(1, x)
    return spherical_y(l - 1, x) - (l + 1) / x * spherical_y(l, x)
end

"""
Analytic phase shift for a finite square well from spherical-Bessel matching.
"""
function analytic_phase_shift(l, V0, a, k)
    E = k^2 / 2
    q = sqrt(2 * (E - V0))
    L_int = q * spherical_j_derivative(l, q * a) / spherical_j(l, q * a)
    num = L_int * spherical_j(l, k * a) - k * spherical_j_derivative(l, k * a)
    den = L_int * spherical_y(l, k * a) - k * spherical_y_derivative(l, k * a)
    return atan(num, den)
end

"""
Scattering length from the product formula. [BT07 Eq.(89)]
"""
function scattering_length(sps)
    return sps.box_size + sum(im / kn for kn in sps.momenta)
end

"""
Distance between two phases, reduced modulo π.
"""
function phase_difference(δ1, δ2)
    return abs(mod(δ1 - δ2 + π / 2, π) - π / 2)
end

@testset "SPS API and eigenvalue count" begin
    a = 1.0
    V = square_well(-112.5, a)

    @test !hasmethod(SPS, Tuple{Integer,Integer,Function,Real,Real})
    @test_throws ArgumentError SPS(8, -1, V, a)

    for l = 0:2
        sps = SPS(12, l, V, a)
        @test length(sps.momenta) == 2 * 12 + l
    end
end

@testset "Phase shift for square-well potential" begin
    a = 1.0
    l = 0
    V0 = -112.5
    V = square_well(V0, a)
    δ_theory = -10.681_639_621_0
    k = 20.0

    function test_phase_shift(N, rtol)
        sps = SPS(N, l, V, a)
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

@testset "Phase shift for nonzero angular momenta" begin
    a = 1.0
    V0 = -112.5
    V = square_well(V0, a)
    k = 20.0

    for l in (1, 2)
        δ_ref = analytic_phase_shift(l, V0, a, k)

        for (N, tol) in ((25, 1e-5), (35, 1e-8))
            sps = SPS(N, l, V, a)
            S = scattering_matrix(sps)(k)
            δ_num = 0.5 * atan(imag(S), real(S))
            @test phase_difference(δ_num, δ_ref) <= tol
        end
    end
end

@testset "Scattering length" begin
    @testset "Rectangular potential" begin
        V = square_well(-112.5, 1.0)
        κ = sqrt(-2 * -112.5)
        α_exact = 1.0 - tan(κ) / κ
        α_paper = 1.057_066_227

        α17 = scattering_length(SPS(17, 0, V, 1.0))
        @test isapprox(real(α17), α_paper; atol = 1e-9, rtol = 0)
        @test abs(imag(α17)) < 1e-10

        α25 = scattering_length(SPS(25, 0, V, 1.0))
        @test isapprox(real(α25), α_exact; atol = 1e-12, rtol = 0)
        @test abs(imag(α25)) < 1e-10
    end

    @testset "Eckart potential (paper Table I)" begin
        eckart(a) = r -> r > a ? 0.0 : -21.0 / cosh(r)^2
        table =
            ((6, 29, 2.458575), (8, 33, 2.450348), (10, 37, 2.450011), (12, 40, 2.450000))

        for (a, N, α_paper) in table
            α = scattering_length(SPS(N, 0, eckart(a), a))
            @test isapprox(real(α), α_paper; atol = 1e-6, rtol = 0)
            @test abs(imag(α)) < 1e-10
        end
    end
end

@testset "Bound-state count for rectangular potential (paper Fig. 1)" begin
    V = square_well(-112.5, 1.0)
    expected = (5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 0)

    for l = 0:11
        sps = SPS(20, l, V, 1.0)
        nbound = count(k -> abs(real(k)) < 1e-8 && imag(k) > 0, sps.momenta)
        @test nbound == expected[l+1]
    end
end

@testset "Free-particle scattering matrix" begin
    a = 1.0
    V = r -> 0.0

    for l = 0:2
        sps = SPS(12, l, V, a)
        S = scattering_matrix(sps)
        for k in (0.5, 1.0, 2.5)
            @test isapprox(S(k), 1.0; atol = 1e-10, rtol = 0)
        end
    end
end
