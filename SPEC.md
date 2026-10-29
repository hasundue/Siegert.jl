# SPEC -- Non-zero angular momentum support

## Goal

`SPS(N, l, V, a)` produces correct one-channel Siegert pseudostates for every
non-negative integer `l`. The scattering matrix and phase shift derived from
those states are correct for `l > 0`, while the existing s-wave results continue
to pass. The obsolete boundary parameter `b` is removed from the public API.

## Background

Current code:

- `src/dvr.jl` already builds a Jacobi DVR basis with `β = 2l`, so the basis
  itself is parameterized by `l`.
- `src/sps.jl` assembles the TON98 matrices and solves a linearized quadratic
  eigenvalue problem. The linearization currently uses `b = -1` and dimension
  `2N`. That is the s-wave boundary condition only.
- `src/spectra.jl` applies the product formula
  `S(k) = exp(-2ika) * prod((kn + k)/(kn - k))` over `sps.momenta`. It will work
  unchanged once `sps.momenta` contains the correct `2N + l` values.

The correct generalization is Batishchev and Tolstikhin, *Phys. Rev. A* **75**,
062704 (2007), "Siegert pseudostate formulation of scattering theory: Nonzero
angular momenta in the one-channel case". The relevant equations are Eq. (20)
(the non-linear eigenvalue problem), Eq. (31) (the exact linearization to
`2N + l`), and Eq. (85) (the product formula for `S_l(k)`).

Decisions already accepted by the user:

- Keep the existing TON98 radial-variable implementation and extend it; do not
  rewrite the package in Batishchev's `u`-variable basis.
- Remove `b` from `SPS`; compatibility with the old five-argument constructor
  is not required.
- `l` remains an integer.
- Keep the existing s-wave tests and add non-zero-`l` tests against an analytic
  square-well reference and a free-particle reference.

## Acceptance criteria

- [ ] `SPS(N, l, V, a)` is the public constructor; the five-argument
      `SPS(N, l, V, a, b)` form is gone.
      verify: `julia --project -e 'using Siegert; V = square_well(-112.5, 1.0); SPS(8, 1, V, 1.0)'`
      exits successfully, and `methods(SPS)` contains no five-argument method.

- [ ] For `l = 0, 1, 2` and `N = 12`, `length(sps.momenta) == 2N + l`.
      verify: `just test sps`

- [ ] Negative `l` is rejected with `ArgumentError`.
      verify: `just test sps`

- [ ] The existing s-wave square-well phase-shift tests still pass after removing
      `b` and still use the same reference value and tolerances.
      verify: `just test sps`

- [ ] For the square well `V0 = -112.5`, `a = 1.0`, `k = 20.0`, the `l = 1` and
      `l = 2` phase shifts agree with the analytic spherical-Bessel matching
      result. At `N = 25` the phase error is at most `1e-5`; at `N = 35` it is at
      most `1e-8`.
      verify: `just test sps`

- [ ] For `V(r) = 0`, `S_l(k) == 1` for `l = 0, 1, 2` and several real `k > 0`
      to numerical tolerance.
      verify: `just test sps`

- [ ] The existing DVR basis tests still pass.
      verify: `just test dvr`

- [ ] Formatting is clean.
      verify: `treefmt -c`

## Out of scope

- Scattering cross sections, physical states, and Green's-function expansions
  from Sec. III of Batishchev and Tolstikhin; only the SPS eigenvalues,
  scattering matrix, and phase shift are in scope.
- Coulomb tails and arbitrary tail corrections. The old `b` parameter was the
  only hook for those; it is being removed.
- Non-integer `l`.
- Numerical guarantees for very large `l`. Zeros of the reverse Bessel
  polynomial are obtained from a companion matrix; Batishchev and Tolstikhin
  report coefficient-based recovery becomes unstable beyond about `l = 28`.
  Do not add a new dependency to extend this range.
- Rewriting the DVR or basis to Batishchev's `u`-variable formulation.
- Fixing stale examples or docstrings unrelated to `SPS` (for example
  `examples/basis.jl` and the `square_well` docstring).
- Multi-channel scattering, non-spherical potentials, and complex scaling.

## Stop conditions

- A required reference cannot be obtained through the publisher's official page
  -> stop and report which reference is missing. Do not substitute an
  unauthorized copy or a secondary summary.
- The same failure twice -> write the cause under `## Blocker`, stop.
- A new dependency is needed -> stop before editing `Project.toml` or
  `Manifest.toml`.
- The s-wave tests cannot be restored to their current tolerances -> stop and
  report the discrepancy.
- The analytic square-well reference and the product formula disagree even as
  `N` grows -> stop; do not weaken the test to hide a formulation error.
- More than three attempts to stabilize the reverse-Bessel zeros -> stop and
  report.
- The work needs to touch anything outside `src/sps.jl`, `src/spectra.jl`,
  `test/test_sps.jl`, and the `SPS` docstring -> stop.

## Decisions

- 2026-10-29 -- Use the Batishchev-Tolstikhin 2007 generalization -- it is the
  paper that explicitly solves the non-zero-`l` one-channel case left open by
  TON98.
- 2026-10-29 -- Keep the current radial variable and DVR, and add the
  `l`-dependent boundary term plus `2N + l` linearization -- smaller diff, and
  the `l = 0` matrix reduces exactly to the current one.
- 2026-10-29 -- Delete `b` -- the user does not need backward compatibility.
- 2026-10-29 -- `l` stays integer -- the formulation and existing API use
  integer angular momentum.
- 2026-10-29 -- Verify against an analytic square well and a free particle --
  both are independent of the SPS implementation.

## Confirmed

- 2026-10-29 -- user: "OK"

## Notes

### Boundary term in the current variables

Current code uses TON98's radial function `R(r)`, while Batishchev uses the
reduced radial function `u(r) = r R(r)`. The outgoing solution in Batishchev's
variables is `u(r) = e_l(kr)`, where `e_l` is defined by the reverse Bessel
polynomial `θ_l(z)`. Its logarithmic derivative at `r = a` is

```text
u'(a) / u(a) = ik - (1/a) * Σ_{p=1..l} z_lp / (ika + z_lp)
```

with `z_lp` the zeros of `θ_l`. Since `u = rR`,

```text
R'(a) / R(a) = u'(a) / u(a) - 1/a
             = ik - 1/a - (1/a) * Σ z_lp / (ika + z_lp).
```

TON98's boundary parameter is defined by `R'(a) = (b/a + ik) R(a)`, so the
correct non-zero-`l` boundary parameter is

```text
b_l(λ) = -1 - Σ_{p=1..l} z_lp / (λa + z_lp),    λ = ik.
```

For `l = 0` the sum is empty and `b = -1`, exactly as today.

### Linearization

Starting from TON98 Eq. (C15),

```text
H̃ c - (b + ika) L c - (k²/2) ρ c = 0,    λ = ik,
```

and substituting `b = b_l(λ)` gives

```text
H̃ c + L c + S L c - λ a L c + (λ²/2) ρ c = 0,
S = Σ z_lp / (λa + z_lp).
```

Introduce `c̃ = λc` and auxiliary scalars `ξ_p`. Let `π = π(a)` so that
`L = π πᵀ`, and choose

```text
ξ_p = - z_lp / (λa + z_lp) * (πᵀ c).
```

Then `(λa + z_lp) ξ_p = -z_lp πᵀc`, hence

```text
λ ξ_p = -(z_lp/a) πᵀc - (z_lp/a) ξ_p,
S L c = -π Σ_p ξ_p.
```

The linearized eigenvalue problem is `M v = λ v` with `v = [c; c̃; ξ]` and

```text
M = [ 0                I                      0        ]
    [ -2ρ⁻¹(H̃ + L)      2aρ⁻¹L                2ρ⁻¹π 1ᵀ ]
    [ -(z_lp/a)πᵀ       0                      -(z_lp/a) I_l ]
```

where the last block is diagonal over `p = 1..l`. The eigenvalue is `λ`; the
SPS momentum is `k = -iλ`. The `c` part of the eigenvector is the first `N`
components.

For `l = 0` this must reduce exactly to the current matrix

```text
[ 0                          I          ]
[ -2ρ⁻¹(H̃ + L)              2aρ⁻¹L     ]
```

which is the current `sps_linearize_qep(..., b = -1)` result.

### Reverse Bessel zeros

Use the recurrence

```text
θ_0(z) = 1
θ_1(z) = z + 1
θ_l(z) = (2l - 1) θ_{l-1}(z) + z² θ_{l-2}(z),   l ≥ 2.
```

The leading coefficient of `θ_l` is 1. Build the monic coefficient vector and
find roots with the eigenvalues of the companion matrix. Return `Float64[]` for
`l = 0`, and `[-1.0]` for `l = 1` (either by the general routine or as a
shortcut). Do not add a polynomial-roots dependency.

### Analytic square-well reference for tests

For `V(r) = V0` for `r < a` and `0` outside, with `E = k²/2` and
`q = sqrt(2(E - V0))`, let `j_l` and `y_l` be the spherical Bessel functions.
The interior logarithmic derivative is

```text
L_int = q * j_l'(qa) / j_l(qa).
```

Matching `u_out(r) ∝ r [cos δ_l j_l(kr) - sin δ_l y_l(kr)]` at `r = a` gives

```text
tan δ_l = (L_int j_l(ka) - k j_l'(ka))
        / (L_int y_l(ka) - k y_l'(ka)).
```

The phase is only defined modulo `π`; compare `S = exp(2iδ)` or reduce the
difference modulo `π`, following the existing s-wave test. `SpecialFunctions`
already provides `besselj` and `bessely`; spherical Bessel functions can be
formed from half-integer orders or from the standard recurrence. No new
dependency is needed.

For the chosen square well (`V0 = -112.5`, `a = 1`, `k = 20`), the reference
phases modulo `π` are approximately

```text
l = 1:  δ ≈ -1.3163168543513333  (mod π)
l = 2:  δ ≈ -1.3106083861911673  (mod π)
```

These values were obtained from the matching formula above; the test should
recompute them rather than hard-code them if practical.

### Files expected to change

- `src/sps.jl`: reverse-Bessel zeros, `l`-dependent linearization, remove `b`,
  update `SPS` docstring.
- `src/spectra.jl`: docstring/comment updates only; the product formula already
  loops over all `sps.momenta`.
- `test/test_sps.jl`: remove the `b` argument, add `2N+l`, analytic `l = 1, 2`,
  free-particle, and negative-`l` tests.
