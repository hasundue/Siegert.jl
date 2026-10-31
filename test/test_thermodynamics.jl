using Test
using Siegert

# DOS electron number used by the thermodynamics neutrality check.
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

@testset "Thermodynamics: identity and DOS neutrality" begin
    Z = 5.0
    aa = average_atom(Z, 10.0, 10.0, 0.2; N = 30, l_max = 6, xc = :ksdt)
    th = thermodynamics(aa)
    T_Ha = aa.temperature_eV / Siegert.HARTREE_TO_EV

    @test aa.converged
    @test all(
        isfinite,
        (
            th.free_energy,
            th.internal_energy,
            th.entropy,
            th.electron_pressure,
            th.ion_ideal_pressure,
            th.total_pressure,
            th.free_energy_electrostatic,
            th.free_energy_xc,
            th.internal_energy_kinetic,
            th.entropy_ks,
            th.entropy_xc,
            th.pressure_xc,
        ),
    )

    # F = U - T_Ha * (S_ks + S_xc)/k_B, with the entropy fields in Ha/K.
    residual =
        th.free_energy - (
            th.internal_energy -
            T_Ha * (th.entropy_ks + th.entropy_xc) / Siegert.KB_HARTREE_PER_K
        )
    @test abs(residual) < 1e-8

    N_dos = dos_electron_number(aa)
    println("thermodynamics: N_dos = ", N_dos, ", relative error = ", (N_dos - Z) / Z)
    @test abs(N_dos - Z) / Z < 1e-6

    # Convenience properties.
    @test th.free_energy_eV ≈ th.free_energy * Siegert.HARTREE_TO_EV
    @test th.internal_energy_eV ≈ th.internal_energy * Siegert.HARTREE_TO_EV
    @test th.entropy_kB ≈ th.entropy / Siegert.KB_HARTREE_PER_K
    @test th.electron_pressure_GPa ≈ th.electron_pressure * Siegert.HARTREE_PER_BOHR3_TO_GPA
    @test th.total_pressure_GPa ≈ th.total_pressure * Siegert.HARTREE_PER_BOHR3_TO_GPA

    # A potential-only XC model must not silently enter the EOS.
    potential_only = average_atom(
        5.0,
        10.0,
        10.0,
        0.2;
        N = 20,
        l_max = 4,
        xc = (n, r) -> zeros(length(n)),
    )
    @test_throws ArgumentError thermodynamics(potential_only)
end

@testset "Thermodynamics: free-electron gas limit" begin
    Z = 1.0
    R = 1.0
    n = Z / (4π * R^3 / 3)
    k_F = (3 * π^2 * n)^(1 / 3)
    ε_F = k_F^2 / 2
    t = 0.01
    T_Ha = t * ε_F
    T_eV = T_Ha * Siegert.HARTREE_TO_EV

    free = fixed_potential_density(
        r -> 0.0,
        Z,
        R,
        T_Ha;
        N = 50,
        l_max = 10,
        density_method = :direct_green,
    )

    # Synthetic uniform free-electron result with V_eff = 0.
    density = fill(n, length(free.grid))
    aa_free = AverageAtomResult(
        density,
        zeros(length(free.grid)),
        free.chemical_potential * Siegert.HARTREE_TO_EV,
        R,
        T_eV,
        0.0,
        0.0,
        0,
        true,
        free.states,
        free.grid,
        free.contour,
        Z,
        NoXC(),
        :direct_green,
        1000,
        2.0,
    )
    th = thermodynamics(aa_free)

    # Sommerfeld expansion to O(t²), sufficient at t = 0.01.
    V = 4π * R^3 / 3
    U_kin = Z * ε_F * (3 / 5) * (1 + (5 * π^2 / 12) * t^2)
    S_kB = Z * (π^2 / 2) * t
    F_kin = U_kin - T_Ha * S_kB
    F_H = 3 * Z^2 / (5 * R)
    # The electron pressure is the continuum/free-electron decomposition:
    # bound/core states do not contribute to the pressure virial.
    P_ref = 2 * U_kin / (3 * V)

    println(
        "free-electron: U_kin = ",
        th.internal_energy_kinetic,
        " vs ",
        U_kin,
        ", S/kB = ",
        th.entropy_kB,
        " vs ",
        S_kB,
        ", F = ",
        th.free_energy,
        " vs ",
        F_H + F_kin,
        ", P = ",
        th.electron_pressure,
        " vs ",
        P_ref,
    )

    @test isapprox(th.internal_energy_kinetic, U_kin; rtol = 1e-4)
    @test isapprox(th.entropy_kB, S_kB; rtol = 1e-3)
    @test isapprox(th.free_energy_electrostatic, F_H; rtol = 1e-6)
    @test isapprox(th.free_energy, F_H + F_kin; rtol = 1e-4)
    @test isapprox(th.internal_energy, F_H + U_kin; rtol = 1e-4)
    @test isapprox(th.electron_pressure, P_ref; rtol = 1e-2)
end

@testset "Thermodynamics: N = 50 vs N = 100" begin
    results = Dict{Int,AverageAtomResult}()
    thermo = Dict{Int,ThermodynamicProperties}()
    for N in (50, 100)
        aa = average_atom(5.0, 10.0, 10.0, 0.2; N = N, l_max = 6, xc = :ksdt)
        @test aa.converged
        results[N] = aa
        thermo[N] = thermodynamics(aa)
    end

    P50 = thermo[50].electron_pressure
    P100 = thermo[100].electron_pressure
    S50 = thermo[50].entropy_kB
    S100 = thermo[100].entropy_kB
    F50 = thermo[50].free_energy
    F100 = thermo[100].free_energy
    U50 = thermo[50].internal_energy
    U100 = thermo[100].internal_energy

    println(
        "thermodynamics N convergence: ΔP/P = ",
        abs(P50 - P100) / abs(P100),
        ", ΔS = ",
        abs(S50 - S100),
        ", ΔF = ",
        abs(F50 - F100),
        ", ΔU = ",
        abs(U50 - U100),
    )
    @test abs(P50 - P100) / abs(P100) < 1e-2
    @test abs(S50 - S100) < 1e-3
    @test abs(F50 - F100) < 1e-4
    @test abs(U50 - U100) < 1e-4
end
