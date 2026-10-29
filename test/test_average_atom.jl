using Test
using Siegert
using LinearAlgebra: I

function eq6_matrix(aa::AverageAtomSPS)
    a = aa.sps.box_size
    l = aa.sps.angular_momentum
    z = Siegert.reverse_bessel_zeros(l)
    boundary = aa.dvr.basis[:, end]
    C = aa.coefficients
    rho = aa.sps.metric
    k = aa.sps.momenta
    M = zeros(ComplexF64, length(k), length(k))
    for i in eachindex(k), j in eachindex(k)
        Pi = a * sum(boundary .* C[:, i])
        Pj = a * sum(boundary .* C[:, j])
        surface =
            1 + sum(
                (zz / ((im * k[i] * a + zz) * (im * k[j] * a + zz)) for zz in z);
                init = zero(ComplexF64),
            )
        M[i, j] =
            (a / 2) * sum(C[:, i] .* (rho * C[:, j])) +
            im * Pi * Pj / (k[i] + k[j]) * surface
    end
    return M
end

function bound_density(states, r, mu, T)
    density = zeros(Float64, length(r))
    for aa in states
        l = aa.sps.angular_momentum
        P = siegert_wavefunction(aa, r)
        for i in eachindex(aa.sps.momenta)
            k = aa.sps.momenta[i]
            if abs(real(k)) < 1e-8 && imag(k) > 0
                epsilon = real(k^2) / 2
                occupation = fermi_dirac(epsilon, mu, T)
                density .+=
                    real.(occupation * 2 * (2l + 1) / (4π) .* (P[:, i] .^ 2 ./ r .^ 2))
            end
        end
    end
    return density
end

relative_l2(a, b) = sqrt(sum((a .- b) .^ 2) / sum(b .^ 2))

# SS23 Eq. (6) / BT07 Eq. (37): the Siegert inner product is bilinear, with
# no complex conjugation. `eq6_matrix` assembles the full normalized Gram
# matrix from the DVR coefficients, metric, and surface term; it must equal I.
@testset "Average atom: Siegert normalization" begin
    V = square_well(-112.5, 1.0)

    # Raw states stay unnormalized for scattering-only use; AA data lives in
    # AverageAtomSPS, so a raw SPS cannot reach the wavefunction layer.
    raw = SPS(12, 0, V, 1.0)
    @test_throws ArgumentError siegert_wavefunction(raw)

    # Include all 2N+l states, hence also the complex resonant pairs.
    for l = 0:2
        aa = AverageAtomSPS(SPS(12, l, V, 1.0))
        @test !isempty(aa.normalizations)
        @test maximum(abs.(eq6_matrix(aa) - I)) < 1e-8
    end
end

# DVR cardinality gives P(r_j) = r_j c_j / kappa_j on the native grid, so this
# checks the FBR/DVR transform and the P = r phi convention. The common grid is
# the l = 0 grid; every partial wave must be evaluable there for Green sums.
@testset "Average atom: wavefunction evaluation" begin
    V = square_well(-112.5, 1.0)
    a = 1.0
    states = [AverageAtomSPS(SPS(12, l, V, a)) for l = 0:2]
    common_r, _ = common_radial_grid(states)

    for aa in states
        P = siegert_wavefunction(aa)
        expected = reshape(aa.sps.grid, :, 1) .* aa.coefficients ./ aa.dvr.coefficients
        @test maximum(abs.(P - expected)) < 1e-10
        @test P[:, 1] == siegert_wavefunction(aa, aa.sps.grid)[:, 1]

        P_common = siegert_wavefunction(aa, common_r)
        @test size(P_common) == (length(common_r), length(aa.sps.momenta))
        @test all(isfinite, P_common)

        # Boundary value used by Eq. (6): P(R) = a pi(1)^T c.
        boundary = a * sum(aa.dvr.basis[:, end] .* aa.coefficients[:, 1])
        @test isapprox(siegert_wavefunction(aa, [a])[1, 1], boundary; rtol = 1e-10)
    end
end

@testset "Average atom: Green's-function residue" begin
    V = square_well(-112.5, 1.0)
    aa = AverageAtomSPS(SPS(20, 0, V, 1.0))
    r = aa.sps.grid
    P = siegert_wavefunction(aa)
    index = findfirst(k -> abs(real(k)) < 1e-8 && imag(k) > 0, aa.sps.momenta)
    @test index !== nothing
    k0 = aa.sps.momenta[index]
    residue = P[:, index] .^ 2 ./ k0
    delta = 1e-8
    G = siegert_green(aa, k0 + delta; r = r)
    @test maximum(abs.(delta .* G .- residue)) / maximum(abs.(residue)) < 1e-6
end

@testset "Average atom: free-electron DOS" begin
    R_values = (1.0, 2.0)
    N_values = (5, 6, 7)
    for N in N_values, R in R_values
        states = [AverageAtomSPS(SPS(N, l, r -> 0.0, R)) for l = 0:5]
        r, weights = common_radial_grid(states)
        for k in (0.5, 1.0, 2.0)
            dos = siegert_dos(states, k; r = r, weights = weights)
            free = (4π * R^3 / 3) * k / π^2
            @test dos > 0
            @test isapprox(dos, free; rtol = 1e-2)
            println("free DOS: N=$N R=$R k=$k dos=$dos free=$free")
        end
    end
end

@testset "Average atom: bound-state limit of the contour" begin
    V = square_well(-112.5, 1.0)
    R = 1.0
    T = 0.1
    mu = -10.0
    states = [AverageAtomSPS(SPS(20, l, V, R)) for l = 0:2]
    r, _ = common_radial_grid(states)
    density = siegert_density(states, mu, T)
    reference = bound_density(states, r, mu, T)
    contour = contour_parameters(mu, T)
    @test contour.E_max == 0.0
    @test contour.N_mat == 0
    @test relative_l2(density, reference) < 1e-10
end

@testset "Average atom: contour invariance" begin
    R = 1.0
    V = r -> r >= R ? 0.0 : -20.0 * (1 - (r / R)^2)^2
    T = 1.0
    mu = 0.5
    states = [AverageAtomSPS(SPS(30, l, V, R)) for l = 0:3]
    reference = siegert_density(states, mu, T)

    more_panels = siegert_density(states, mu, T; n_vert = 2, n_horiz = 40, n_return = 40)
    other_EI = siegert_density(states, mu, T; E_I = 1.2)
    lower_E_min = siegert_density(states, mu, T; E_min = -0.5)

    @test relative_l2(more_panels, reference) < 1e-6
    @test relative_l2(other_EI, reference) < 1e-6
    @test relative_l2(lower_E_min, reference) < 1e-6
end

@testset "Average atom: neutrality and fixed potential" begin
    V = square_well(-5.0, 1.0)
    Z = 5.0
    R = 1.0
    T = 0.2
    result = fixed_potential_density(V, Z, R, T; N = 30, l_max = 3)

    @test isfinite(result.chemical_potential)
    @test abs(result.neutrality_residual) / Z < 1e-8
    @test minimum(result.density) > -1e-12
    @test result.contour.N_mat >= 1
    @test result.contour.E_max > 0

    r = result.grid
    weights = (R / 2) .* result.states[1].dvr.weights
    electron_number = 4π * sum(volume_weights(r, weights) .* result.density)
    @test isapprox(electron_number, Z; rtol = 1e-8)

    nbound = count(
        k -> abs(real(k)) < 1e-8 && imag(k) > 0,
        reduce(vcat, (aa.sps.momenta for aa in result.states)),
    )
    @test nbound > 0
end

@testset "Average atom: basis convergence" begin
    V = square_well(-5.0, 1.0)
    Z = 5.0
    R = 1.0
    T = 0.3
    l_max = 3
    results = Dict{Int,AADensity}()
    for N in (50, 100, 150)
        result = fixed_potential_density(V, Z, R, T; N = N, l_max = l_max)
        results[N] = result
        println("convergence: N=$N mu=$(result.chemical_potential)")
    end

    mus = [results[N].chemical_potential for N in (50, 100, 150)]
    @test maximum(mus) - minimum(mus) < 1e-4

    r_ref = collect(range(1e-3, R; length = 100))
    densities = [
        siegert_density(
            results[N].states,
            results[N].chemical_potential,
            results[N].temperature;
            r = r_ref,
            weights = ones(length(r_ref)),
        ) for N in (50, 100, 150)
    ]
    for i = 1:2
        change = relative_l2(densities[i+1], densities[i])
        println("convergence: density relative L2 change $(i) = $change")
        @test change < 1e-2
    end
end
