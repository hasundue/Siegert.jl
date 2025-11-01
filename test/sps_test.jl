using Test
using Siegert
using LinearAlgebra: norm
using Printf: @sprintf
using Statistics: median, quantile
using Logging: @debug, @warn

@testset "SPS eigenvalue tests (square well)" begin
    a = 1.0
    l = 0
    V0 = -112.5
    V = _ -> V0
    # Paper baseline: 2N = 38 → N = 19
    N = 19
    ks, C, meta = sps_solve(N, l, a, V; b = 0.0)

    # Filter those close to the imaginary axis and (optionally) log for Im>0 and Im<0
    axis_eps = 1e-8
    ks_axis = filter(k -> abs(real(k)) < axis_eps, ks)
    ks_pos = sort(filter(k -> imag(k) > 0, ks_axis), by = k -> imag(k), rev = true)
    ks_neg = sort(filter(k -> imag(k) < 0, ks_axis), by = k -> imag(k))

    @debug "Near-imaginary-axis eigenvalues (|Re k| < $(axis_eps))"
    @debug "Im > 0 (largest → smallest)" ks_pos=ks_pos
    @debug "Im < 0 (most negative → least negative)" ks_neg=ks_neg

    # Basic: we should have 2N eigenvalues
    @test length(ks) == 2N

    # Sanity 1: existence of near-imaginary eigenvalues (Re k ≈ 0)
    num_near_imag = count(abs(real(k)) < 1e-8 for k in ks)
    @test num_near_imag ≥ 2

    # Symmetry pairing metric: k ↔ −conj(k)
    eps = 1e-10
    right = [k for k in ks if real(k) > eps]
    if !isempty(right)
        dists = [minimum(abs(kj + conj(kr)) for kj in ks) for kr in right]
        med = median(dists)
        p95 = quantile(dists, 0.95)
        @debug "Pairing distances |k + conj(k′)|" median=med p95=p95
        @test med < 1e-8
        @test p95 < 1e-6
    else
        @warn "No right-half-plane eigenvalues to pair; skipping pairing checks"
    end
end
