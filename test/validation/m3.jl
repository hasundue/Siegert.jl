# M3c validation: SS23 figures 1-6 and the low-temperature Al/Na/Fe cases.
#
# Usage:
#   julia --project test/validation/m3.jl          # SS23 Figs. 1-6 subset
#   julia --project test/validation/m3.jl lowt     # low-temperature cases
#   julia --project test/validation/m3.jl all      # both
#
# The pressure in SS23 Fig. 6 is the Tartarus continuum/free-electron pressure
# (bound states excluded from the pressure virial). The reference values below
# were digitized from the vector curves in the SS23 source `pressure.pdf`.

using Test
using Siegert
using Serialization

const N_VALIDATION = 50
const LMAX_AL_FE = 8
const LMAX_C = 6
const LMAX_LOWT = 10

# Digitized SS23 Fig. 6: (T eV, rho g/cm^3) => pressure in Mbar.
const SS23_FIG6 = Dict(
    (10.0, 1e-3) => 0.00229804,
    (10.0, 1.0) => 1.14748,
    (10.0, 100.0) => 3432.52,
    (100.0, 1e-3) => 0.0474914,
    (100.0, 1.0) => 33.7348,
    (100.0, 100.0) => 5054.06,
    (1000.0, 1e-3) => 0.481545,
    (1000.0, 1.0) => 473.747,
    (1000.0, 100.0) => 45095.0,
)

const AA_CACHE = Dict{NTuple{6,Float64},AverageAtomResult}()
const CACHE_PATH =
    normpath(joinpath(@__DIR__, "..", "..", "scratch", "m3_validation_cache.jls"))

if isfile(CACHE_PATH)
    try
        loaded = deserialize(CACHE_PATH)
        loaded isa Dict{NTuple{6,Float64},AverageAtomResult} && merge!(AA_CACHE, loaded)
        println("loaded $(length(AA_CACHE)) cached average-atom results")
    catch err
        @warn "could not load validation cache" exception = err
    end
end

function get_aa(
    Z::Real,
    A::Real,
    rho::Real,
    T::Real;
    N::Integer = N_VALIDATION,
    l_max::Integer = LMAX_AL_FE,
)
    key = (Float64(Z), Float64(A), Float64(rho), Float64(T), Float64(N), Float64(l_max))
    if !haskey(AA_CACHE, key)
        println("running: Z=$Z A=$A rho=$rho T=$T eV, N=$N, l_max=$l_max ...")
        t0 = time()
        aa = average_atom(Z, A, rho, T; N = N, l_max = l_max, xc = :ksdt)
        println("  done in $(round(time() - t0; digits = 1)) s, converged=$(aa.converged)")
        AA_CACHE[key] = aa
        try
            serialize(CACHE_PATH, AA_CACHE)
        catch err
            @warn "could not write validation cache" exception = err
        end
    end
    return AA_CACHE[key]
end

function partial_dos(aa::AverageAtomResult, k::Real, l::Integer)
    return siegert_dos(aa.states[l+1], k; r = aa.grid, weights = Siegert._dr_weights(aa))
end

function total_dos(aa::AverageAtomResult, k::Real)
    return siegert_dos(aa.states, k; r = aa.grid, weights = Siegert._dr_weights(aa))
end

function bound_momenta(aa::AverageAtomResult, l::Integer; tolerance = 1e-8)
    ks = aa.states[l+1].sps.momenta
    return [k for k in ks if Siegert._is_bound_state(k, tolerance)]
end

function resonance_momenta(aa::AverageAtomResult, l::Integer; tolerance = 1e-8)
    ks = aa.states[l+1].sps.momenta
    return [k for k in ks if real(k) > tolerance && imag(k) < -tolerance]
end

function select_resonance(aa::AverageAtomResult, l::Integer, lo::Real, hi::Real)
    resonances = filter(k -> lo < real(k) < hi, resonance_momenta(aa, l))
    isempty(resonances) && error("no resonance with k in [$lo, $hi] for l=$l")
    return resonances[argmin(abs.(imag.(resonances)))]
end

function dos_peak(f, lo::Real, hi::Real; step::Real = 0.005)
    best_k = lo
    best = -Inf
    k = lo
    while k <= hi
        value = f(k)
        if value > best
            best = value
            best_k = k
        end
        k += step
    end
    return best_k, best
end

function dos_electron_number(aa::AverageAtomResult; dos_cutoff = 1e-10)
    T = aa.temperature_eV / Siegert.HARTREE_TO_EV
    μ = aa.chemical_potential_eV / Siegert.HARTREE_TO_EV
    bound = Siegert._bound_dos(aa)
    N = sum(weight * fermi_dirac(ε, μ, T) for (ε, weight) in bound)
    R = aa.radius_bohr
    free_setup = (aa.grid, Siegert._dr_weights(aa), R)
    direct_setup =
        aa.density_method === :direct_green ? Siegert._direct_dos_setup(aa) : nothing
    l_max = length(aa.states) - 1
    E_high = μ + T * log(1 / dos_cutoff - 1)
    if E_high > aa.contour.E_min
        points, measures = Siegert._energy_quadrature(aa.contour.E_min, E_high, μ, T)
        for (ε, w) in zip(points, measures)
            χ = Siegert._continuum_dos(aa, ε, direct_setup, free_setup, l_max)
            N += w * fermi_dirac(ε, μ, T) * χ
        end
    end
    return N
end

function all_state_virial(aa::AverageAtomResult, th::ThermodynamicProperties)
    V = 4π * aa.radius_bohr^3 / 3
    return (2 * th.internal_energy_kinetic + th.free_energy_electrostatic) / (3V) +
           th.pressure_xc
end

function validate_ss23()
    @testset "SS23 Fig. 1: Al l=0 eigenvalues and DOS" begin
        dos10 = nothing
        dos20 = nothing
        dos30 = nothing
        for T in (10.0, 20.0, 30.0)
            aa = get_aa(13.0, 27.0, 2.7, T; l_max = LMAX_AL_FE)
            @test aa.converged
            bound = bound_momenta(aa, 0)
            bound_im = sort(imag.(bound); rev = true)

            # The 2s bound state is present for all three temperatures.
            @test any(k -> 2.0 < imag(k) < 4.0, bound)
            if T == 30.0
                # The 3s state is bound at 30 eV.
                @test any(k -> 0.05 < imag(k) < 0.8, bound)
            else
                # At 10/20 eV the 3s state is anti-bound, just below the real axis.
                ks = aa.states[1].sps.momenta
                @test any(k -> abs(real(k)) < 1e-8 && -0.8 < imag(k) < -0.01, ks)
            end

            dos = k -> partial_dos(aa, k, 0)
            if T == 10.0
                dos10 = dos
            elseif T == 20.0
                dos20 = dos
            else
                dos30 = dos
            end
            peak_k, peak = dos_peak(dos, 0.02, 1.5; step = 0.01)
            println(
                "Fig. 1: Al T=$T eV, bound l=0 Im(k) = ",
                bound_im,
                ", l=0 DOS peak at k=",
                peak_k,
                " (",
                peak,
                ")",
            )
            @test peak > 3.0
        end
        # The 3s anti-bound state is closest to the axis at 20 eV and gives the
        # largest low-k DOS enhancement.
        @test dos20(0.1) > dos10(0.1)
        @test dos20(0.1) > dos30(0.1)
    end

    @testset "SS23 Fig. 2: Fe l=2 resonance and DOS" begin
        aa = get_aa(26.0, 56.0, 8.0, 1.0; l_max = LMAX_AL_FE)
        @test aa.converged
        k_res = select_resonance(aa, 2, 0.5, 1.2)
        peak_k, peak = dos_peak(k -> partial_dos(aa, k, 2), 0.5, 1.2; step = 0.002)
        println(
            "Fig. 2: Fe l=2 nearest resonance k=",
            k_res,
            ", l=2 DOS peak at k=",
            peak_k,
            " (",
            peak,
            ")",
        )
        @test 0.7 < real(k_res) < 1.0
        @test abs(imag(k_res)) < 0.1
        @test 0.7 < peak_k < 1.0
        @test peak > 40.0
    end

    @testset "SS23 Fig. 3: Fe effective potential" begin
        aa = get_aa(26.0, 56.0, 8.0, 1.0; l_max = LMAX_AL_FE)
        r = aa.grid
        V = aa.potential
        V_l2 = V .+ 1.0 ./ (r .^ 2)
        inner = (r .> 0.6) .& (r .< 3.0)
        V_min = minimum(V[inner])
        well = findmin(V_l2[inner])[1]
        barrier = maximum(V_l2[(r .> 1.0) .& (r .< 4.0)])
        k_res = select_resonance(aa, 2, 0.5, 1.2)
        ε_res = real(k_res^2) / 2
        println(
            "Fig. 3: Fe V_eff min=",
            V_min,
            " Ha, l=2 well=",
            well,
            " Ha, barrier=",
            barrier,
            " Ha, resonance ε=",
            ε_res,
            " Ha",
        )
        @test V_min < -3.0
        @test well < -2.0
        @test barrier > 0.1
        @test 0.2 < ε_res < 0.5
        @test abs(V[end]) < 1e-8
    end

    @testset "SS23 Figs. 4/5: C eigenvalue tracks" begin
        one_s_im = Float64[]
        p_bound_low = 0
        p_res_high = 0
        for rho in (1e-3, 1.0, 100.0)
            aa = get_aa(6.0, 12.0, rho, 10.0; l_max = LMAX_C)
            @test aa.converged
            bound0 = bound_momenta(aa, 0)
            bound1 = bound_momenta(aa, 1)
            res1 = resonance_momenta(aa, 1)
            if !isempty(bound0)
                push!(one_s_im, maximum(imag.(bound0)))
            end
            rho == 1e-3 && (p_bound_low = length(bound1))
            rho == 100.0 && (p_res_high = length(res1))
            println(
                "Figs. 4/5: C rho=$rho, l=0 bound Im(k)=",
                sort(imag.(bound0); rev = true),
                ", l=1 bound count=",
                length(bound1),
                ", l=1 resonance count=",
                length(res1),
            )
        end
        # The 1s state remains bound and moves toward the real axis as density
        # increases; at low density several p states are bound, while at high
        # density pressure-ionized p states appear as resonances.
        @test length(one_s_im) >= 2
        @test issorted(one_s_im; rev = true)
        @test p_bound_low >= 2
        @test p_res_high >= 1
    end

    @testset "SS23 Fig. 6: C excess pressure" begin
        for T in (10.0, 100.0, 1000.0), rho in (1e-3, 1.0, 100.0)
            aa = get_aa(6.0, 12.0, rho, T; l_max = LMAX_C)
            @test aa.converged
            th = thermodynamics(aa)
            P = th.electron_pressure_GPa / 100
            P_vir = all_state_virial(aa, th) * Siegert.HARTREE_PER_BOHR3_TO_GPA / 100
            P_ref = SS23_FIG6[(T, rho)]
            log_error = abs(log10(P / P_ref))
            println(
                "Fig. 6: C rho=$rho T=$T eV: P=",
                round(P; sigdigits = 5),
                " Mbar, ref=",
                P_ref,
                " Mbar, rel=",
                round((P - P_ref) / P_ref; sigdigits = 3),
                ", all-state virial=",
                round(P_vir; sigdigits = 5),
                " Mbar",
            )
            @test isfinite(P)
            # Log-scale figure: 0.1 dex is well inside the line width/resolution.
            @test log_error < 0.1
        end
    end
end

function validate_lowt()
    @testset "Low-temperature Al/Na/Fe" begin
        cases = (
            (13.0, 27.0, 2.7, 0.1, "Al"),
            (11.0, 23.0, 0.97, 0.1, "Na"),
            (26.0, 56.0, 7.0, 0.2, "Fe"),
        )
        for (Z, A, rho, T, name) in cases
            aa = get_aa(Z, A, rho, T; l_max = LMAX_LOWT)
            th = thermodynamics(aa)
            N_dos = dos_electron_number(aa)
            println(
                "lowT $name: conv=",
                aa.converged,
                ", mu=",
                aa.chemical_potential_eV,
                " eV, N_dos=",
                N_dos,
                ", rel=",
                (N_dos - Z) / Z,
                ", min(n)=",
                minimum(aa.density),
                ", P=",
                th.electron_pressure_GPa / 100,
                " Mbar",
            )
            @test aa.converged
            @test minimum(aa.density) > -1e-8
            # The strict 1e-6 criterion applies to the SS23 reference cases. At
            # low T, near-threshold anti-bound states make the real-axis DOS
            # neutrality residual ~1e-5 (still below the 1e-4 bonus criterion).
            @test abs(N_dos - Z) / Z < 1e-4
            @test all(
                isfinite,
                (
                    th.free_energy,
                    th.internal_energy,
                    th.entropy,
                    th.electron_pressure,
                    th.total_pressure,
                ),
            )
        end
    end
end

if isempty(ARGS)
    validate_ss23()
elseif ARGS == ["lowt"]
    validate_lowt()
elseif ARGS == ["all"]
    validate_ss23()
    validate_lowt()
else
    error("usage: julia --project test/validation/m3.jl [lowt|all]")
end
