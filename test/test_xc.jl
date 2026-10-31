using Test
using Siegert
using ForwardDiff

# Reference values from libxc 6.2.2, functional LDA_XC_CORRKSDT (unpolarized),
# evaluated with `xc_lda_exc_vxc` at the given (r_s, t). These are the
# independent check of the closed-form implementation.
const LIBXC_CORRKSDT = [
    (0.1, 1e-6, -4.6916110446916024, -6.2271105812873095),
    (0.5, 1e-4, -0.98981607105169678, -1.3024280315469945),
    (1.0, 0.01, -0.51694581855831767, -0.6766119438149989),
    (1.0, 0.1, -0.51630297469999575, -0.67652141060258109),
    (1.0, 1.0, -0.39545833525768309, -0.58033831616220022),
    (1.0, 10.0, -0.1450883843575583, -0.22315379498711599),
    (2.0, 0.1, -0.27419719604275111, -0.35681022460386347),
    (2.0, 1.0, -0.228137840887953, -0.32497560715824225),
    (4.0, 0.5, -0.14249389684781608, -0.19052840977016833),
    (10.0, 0.1, -0.064652555891770033, -0.083911571127013557),
    (0.1, 10.0, -0.61090090510284678, -1.0092300453311944),
    (10.0, 10.0, -0.035629690935068757, -0.051768282484996569),
    (3.0, 2.0, -0.13879540195063295, -0.20007553109846865),
    (0.5, 0.5, -0.84705359224311394, -1.2198679524245333),
    (5.0, 0.05, -0.12015308565669831, -0.15590502514326382),
    (0.8, 3.0, -0.30363735198665731, -0.46730402018662753),
    (1.5, 0.001, -0.35589485306675545, -0.46443171578604481),
    (0.2, 0.2, -2.2617006930603241, -3.1260201996369554),
]

function density_temperature(rs, t)
    n = 3 / (4π * rs^3)
    t_f = (3 * π^2 * n)^(2 / 3) / 2
    return n, t * t_f
end

relative_l2(a, b) = sqrt(sum((a .- b) .^ 2) / sum(b .^ 2))

@testset "XC: corrKSDT against libxc and limits" begin
    for (rs, t, f_ref, v_ref) in LIBXC_CORRKSDT
        n, T = density_temperature(rs, t)
        f = only(xc_free_energy(KSDT(), [n], T))
        v = only(xc_potential(KSDT(), [n], [0.0], T))
        @test isapprox(f, f_ref; rtol = 1e-12, atol = 1e-14)
        @test isapprox(v, v_ref; rtol = 1e-12, atol = 1e-14)
    end

    # t -> 0: the tanh(1/t) factors reach 1 and exp(-c3/t) vanishes.
    for rs in (0.1, 0.5, 1.0, 3.0, 10.0)
        f0, _, _ = Siegert._corrksdt_f_rs_ft(rs, 0.0)
        f_small, _, _ = Siegert._corrksdt_f_rs_ft(rs, 1e-12)
        @test isfinite(f0)
        @test isapprox(f_small, f0; rtol = 1e-11, atol = 1e-14)
    end

    # High-T leading asymptote from the fit:
    # sqrt(t) f -> -b3 / (b5 sqrt(rs)).
    for rs in (0.1, 1.0, 10.0)
        f, _, _ = Siegert._corrksdt_f_rs_ft(rs, 1e10)
        expected = -Siegert._KSDT_B[3] / (Siegert._KSDT_B[5] * sqrt(rs))
        @test isapprox(sqrt(1e10) * f, expected; rtol = 1e-4)
    end

    @test xc_potential(nothing, [1.0], [0.0], 0.1) == [0.0]
    @test xc_free_energy(nothing, [1.0], 0.1) == [0.0]
end

@testset "XC: analytic derivatives against ForwardDiff" begin
    max_rs_error = 0.0
    max_t_error = 0.0
    for rs in 10 .^ range(-1, 1; length = 16)
        for t in 10 .^ range(-6, 1; length = 16)
            _, f_rs, f_t = Siegert._corrksdt_f_rs_ft(rs, t)
            f_rs_fd =
                ForwardDiff.derivative(x -> first(Siegert._corrksdt_f_rs_ft(x, t)), rs)
            f_t_fd = ForwardDiff.derivative(x -> first(Siegert._corrksdt_f_rs_ft(rs, x)), t)
            rs_error = abs(f_rs - f_rs_fd) / max(abs(f_rs_fd), 1e-300)
            t_error = abs(f_t - f_t_fd) / max(abs(f_t_fd), 1e-300)
            max_rs_error = max(max_rs_error, rs_error)
            max_t_error = max(max_t_error, t_error)
        end
    end
    println(
        "XC derivatives: max relative error f_rs = ",
        max_rs_error,
        ", f_t = ",
        max_t_error,
    )
    @test max_rs_error < 1e-6
    @test max_t_error < 1e-6

    # Approaching t -> 0, including the derivative.
    for rs in (0.1, 1.0, 10.0), t in (1e-6, 1e-8)
        _, _, f_t = Siegert._corrksdt_f_rs_ft(rs, t)
        f_t_fd = ForwardDiff.derivative(x -> first(Siegert._corrksdt_f_rs_ft(rs, x)), t)
        @test abs(f_t - f_t_fd) / max(abs(f_t_fd), 1e-300) < 1e-6
    end
end

@testset "XC: V_xc = f + n df/dn" begin
    max_error = 0.0
    for rs in (0.1, 0.5, 1.0, 3.0, 10.0), t in (1e-6, 1e-4, 0.01, 0.1, 1.0, 10.0)
        n, T = density_temperature(rs, t)
        f_of_n = nn -> only(xc_free_energy(KSDT(), [nn], T))
        v_fd = f_of_n(n) + n * ForwardDiff.derivative(f_of_n, n)
        v = only(xc_potential(KSDT(), [n], [0.0], T))
        max_error = max(max_error, abs(v - v_fd) / max(abs(v), 1e-300))
    end
    println("XC potential: max relative error vs f + n df/dn = ", max_error)
    @test max_error < 1e-6
end

@testset "XC: average_atom with corrKSDT" begin
    result = average_atom(5.0, 10.0, 10.0, 0.2; N = 40, l_max = 4, xc = :ksdt)
    @test result.converged
    @test result.self_consistency_residual < 1e-4
    @test abs(result.neutrality_residual) / 5.0 < 1e-8
    @test abs(result.potential[end]) < 1e-8
    @test result.xc isa KSDT
    @test result.nuclear_charge == 5.0

    # `xc = nothing` stays the M2 zero-XC result, and the NoXC object is
    # equivalent to it.
    base = average_atom(5.0, 10.0, 10.0, 0.2; N = 40, l_max = 4, xc = nothing)
    noxc = average_atom(5.0, 10.0, 10.0, 0.2; N = 40, l_max = 4, xc = NoXC())
    @test relative_l2(noxc.density, base.density) < 1e-14
    @test maximum(abs.(noxc.potential .- base.potential)) < 1e-14

    sugar = average_atom(5.0, 10.0, 10.0, 0.2; N = 40, l_max = 4, xc = KSDT())
    @test relative_l2(sugar.density, result.density) < 1e-14
    @test maximum(abs.(sugar.potential .- result.potential)) < 1e-14

    # A raw callable is potential-only and has no free energy.
    callable(n, r) = zeros(length(n))
    wrapped = average_atom(5.0, 10.0, 10.0, 0.2; N = 40, l_max = 4, xc = callable)
    @test wrapped.xc isa Siegert.PotentialOnlyXC
    @test xc_free_energy(wrapped.xc, wrapped.density, 0.01) === nothing
end
