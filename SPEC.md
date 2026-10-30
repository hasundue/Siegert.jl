# SPEC -- Self-consistent average atom (M2)

## Goal

M2 turns the M1 fixed-potential density layer into a self-consistent
finite-temperature average atom. Given `Z`, mass number `A`, mass density
`rho_m`, and temperature `T`, it returns the self-consistent electron density,
effective potential, chemical potential, neutrality residual, and convergence
information.

The loop uses the Hartree potential and a pluggable exchange-correlation
potential. M2 ships with `V_xc = 0` by default. Free-electron subtraction is a
separate M2 sub-increment. KSDT exchange-correlation, pressure, entropy, and
the equation of state remain M3.

## Background

- `src/average_atom.jl` contains the M1 layer: `NormalizedSPS`, the Siegert
  wavefunctions and Green's function, `siegert_density`, `chemical_potential`,
  `fixed_potential_density`, `Contour`, and `AADensity`.
- `fixed_potential_density(V, Z, R, T; N, l_max)` computes `n_e(r)` and `μ` for
  a fixed `V`. It is the input to M2.
- SS23 Eqs. (10)-(13) define the finite-temperature density, the neutrality
  condition, and the contour + Matsubara density. SS23 Eq. (11) defines
  `V_eff` from `-Z/r`, the Hartree potential, and
  `V_xc[n(r)] - V_xc[n(R)]`.
- S15 Appendix B defines the contour construction. S15 Eq. (13)/(14) is the
  free-electron subtraction reference for the angular-momentum sum.
- S15 explicitly uses the KSDT finite-temperature exchange-correlation free
  energy. SS23 writes only `V_xc` and finite-temperature DFT, but it belongs to
  the same S15/Tartarus lineage. M2 does not implement KSDT.
- Internal units are Hartree and Bohr. The public API accepts `T` in eV and
  `rho_m` in g/cm^3, and returns `chemical_potential_eV` in eV.

## Scope

### M2a: self-consistent Hartree loop

Implement a fixed-point loop over `V_eff`:

1. Start from `V_eff(r) = -Z/r + Z/R` for `0 < r < R`, and `V_eff(r) = 0` for
   `r >= R`. This is the simplest initial guess that satisfies `V_eff(R) = 0`.
2. Compute `n_e(r)` and `μ` with the M1 density layer for the current `V_eff`.
3. Construct the Hartree potential
   `V_H(r) = 4π[1/r ∫_0^r r'^2 n_e dr' + ∫_r^R r' n_e dr']`.
4. Construct
   `V_eff,new(r) = -Z/r + V_H(r) + V_xc[n_e(r)] - V_xc[n_e(R)]`.
5. Mix:
   `V_eff,next = (1 - mixing) V_eff + mixing V_eff,new`.
6. Repeat until convergence or `max_iterations`.

Default `mixing = 0.4`, `max_iterations = 100`. No Anderson or other
acceleration in M2.

Convergence requires

- `max_r |V_eff,next(r) - V_eff(r)| < 1e-4` Ha, and
- `|4π ∫_0^R r^2 n_e dr - Z| / Z < 1e-8`.

The self-consistency threshold is `1e-4` Ha because the M1 contour + Matsubara
density carries a Float64 noise floor of ~1e-5 relative to `max |n_e|`; see
`## Blocker`. The boundary and neutrality thresholds stay at `1e-8`.

If convergence is not reached, return a result with `converged = false` and
emit a warning. The result must contain the final density, potential,
residuals, and iteration count. Tests use cases that converge.

### M2b: pluggable exchange-correlation

`average_atom` accepts a keyword `xc = nothing`. `nothing` means zero
exchange-correlation. Otherwise `xc` is a callable

```julia
xc(n::AbstractVector, r::AbstractVector) -> AbstractVector
```

returning `V_xc[n](r)` on the grid. The boundary subtraction
`V_xc[n(R)]` is applied by the driver, so the supplied function must return
the unsubtracted `V_xc[n](r)`.

The interface is required and is tested with an analytic test functional.
A production simple XC (for example finite-temperature exchange-only) is
optional and is not required for M2 completion. KSDT is M3.

### M2c: free-electron subtraction

Implement free-electron subtraction as an independent sub-increment, enabled
by a keyword such as `free_electron_subtraction = true`.

The intent is the S15 Eq. (13)/(14) construction: subtract the free-electron
Green function from the Siegert Green function in the angular-momentum sum and
add back the free-electron contribution, so that the total density is unchanged
while the `l` sum converges faster. SS23 uses the related device of treating
`l > l_max` as free electrons.

M2c must be validated first at fixed potential using the M1 density layer, then
enabled in the self-consistent loop. If the free reference cannot be defined
consistently with the finite cell, stop and report; do not fit a correction.

## API

New public names:

- `average_atom(Z, A, rho_m, T; N=50, l_max=10, xc=nothing, mixing=0.4,
  max_iterations=100, free_electron_subtraction=false, ...) ->
  AverageAtomResult`
- `AverageAtomResult`: fields
  - `density::Vector{Float64}` in bohr^-3
  - `potential::Vector{Float64}` in Ha, the final `V_eff`
  - `chemical_potential_eV::Float64`
  - `radius_bohr::Float64`
  - `temperature_eV::Float64`
  - `neutrality_residual::Float64` (`4π ∫ r² n_e dr - Z`)
  - `self_consistency_residual::Float64` (`max |ΔV_eff|`)
  - `iterations::Int`
  - `converged::Bool`
  - `states::Vector{NormalizedSPS}`
  - `grid::Vector{Float64}` in bohr
  - `contour::Contour`

`fixed_potential_density` and the M1 API remain available.

## Units and conversions

Public inputs:

- `T` in eV
- `rho_m` in g/cm^3
- `A` in atomic mass units (u)
- `Z` as a positive number

Internal conversions:

- `T_Ha = T_eV / 27.211386245988`
- `R_m = (3 A m_u / (4π rho_m))^(1/3)`
- `R_Bohr = R_m / a_0`
- `chemical_potential_eV = chemical_potential_Ha * 27.211386245988`

Constants: `m_u = 1.66053906660e-27` kg, `a_0 = 5.29177210903e-11` m,
`1 g/cm^3 = 1000 kg/m^3`. These are defined in one place and tested by a unit
conversion test.

## Acceptance criteria

- [ ] With `xc = nothing`, a test case (for example `Z = 5`, `A = 10`,
      `rho_m = 10`, `T = 0.2` eV) converges with `converged = true`,
      `self_consistency_residual < 1e-4`, and
      `abs(neutrality_residual) / Z < 1e-8`.
      verify: `just test average_atom`

- [ ] At convergence, `V_eff(R) = 0` within `1e-8` Ha and the density is
      non-negative within roundoff.
      verify: `just test average_atom`

- [ ] A pluggable analytic test `V_xc` runs through the loop, changes
      `V_eff`, and still converges; the supplied `V_xc` is not subtracted twice.
      verify: `just test average_atom`

- [ ] `average_atom(Z, A, rho_m, T)` converts eV, g/cm^3, and u to Hartree and
      Bohr correctly; the same physical case called with the equivalent
      `R_Bohr` and `T_Ha` gives the same density and `μ`.
      verify: `just test average_atom`

- [ ] `N = 30, 40, 50` on a self-consistent test case gives stable `μ` and
      density. Report the actual numbers in the test output. Thresholds:
      `max |Δμ| < 1e-2` eV and density relative `L2` change `< 1e-2`. The
      `μ` threshold is set by the M1 density noise floor (~1e-4 Ha = 2.7e-3
      eV) and the physical target of ~0.01% in `μ`; it is not a 1e-8 target.
      The fixed-potential `N = 50, 100, 150` sweep remains covered by the M1
      basis-convergence test.
      verify: `just test average_atom`

- [ ] Free-electron subtraction at fixed potential agrees with the
      no-subtraction density at a larger `l_max` to `1e-6` relative `L2`, and
      its on/off results agree at the same large `l_max`.
      verify: `just test average_atom`

- [ ] Existing M1 and DVR/SPS tests still pass.
      verify: `just test dvr sps` or `just test`

- [ ] Formatting is clean.
      verify: `treefmt -c`

## Out of scope

- KSDT exchange-correlation. KSDT is M3.
- Pressure, entropy, free energy, equation of state, and SS23/Tartarus
  comparisons. These are M3.
- A required production simple XC. The pluggable interface is required; a
  simple production model is optional.
- Relativistic effects, non-spherical potentials, multi-channel scattering,
  and non-integer `l`.
- Off-diagonal Green's function `G_l(r,r';k)` for `r != r'`.
- Anderson mixing, Pulay mixing, and other convergence accelerators.
- Changing the M1 public API semantics. `fixed_potential_density` remains.
- Thomas-Fermi or other elaborate initial guesses; M2 uses `-Z/r + Z/R`.

## Stop conditions

- The same failure happens twice -> write the cause under `## Blocker`, stop.
- The self-consistent loop does not converge for the test case with mixing in
  `0.3-0.5` and `max_iterations = 100` -> stop and report; do not add a new
  accelerator or fit a mixing schedule without agreement.
- The free-electron subtraction reference cannot be defined consistently with
  the finite cell -> stop and report; do not fit a correction factor.
- A unit conversion or boundary-value convention is ambiguous -> stop and
  report.
- A new dependency is needed but its source/API cannot be checked, or it is
  not a well-known/de facto package -> stop and ask.
- The work needs changes outside `src/average_atom.jl`, `src/Siegert.jl`,
  `src/sps.jl`, `src/dvr.jl`, `test/test_average_atom.jl`, `README.md`, and
  related docstrings -> stop.

## Decisions

- 2026-10-29 -- M2 is split into M2a (self-consistent Hartree with the
  required pluggable `V_xc` interface), M2b (optional production simple XC),
  and M2c (free-electron subtraction as an independent sub-increment).
  KSDT/pressure/EOS are M3.
- 2026-10-29 -- `V_xc` is a pluggable potential-level callable; M2 default is
  zero. The pluggable interface is required, but a production simple XC is
  optional.
- 2026-10-29 -- Public `average_atom` accepts `T` in eV and `rho_m` in g/cm^3;
  `chemical_potential_eV` is returned in eV. Density is bohr^-3 and radius is
  Bohr.
- 2026-10-29 -- Default mixing is 0.4, maximum iterations 100, convergence
  `max |ΔV_eff| < 1e-8` Ha and relative neutrality `< 1e-8`. No Anderson
  mixing in M2.
- 2026-10-30 -- The self-consistency threshold is relaxed from `1e-8` to
  `1e-4` Ha (decision A). The M1 density layer's Float64 noise floor is
  ~1e-5, so `1e-8` is not attainable without changing the M1 density
  representation. The boundary `V_eff(R) = 0` within `1e-8` Ha and the
  relative neutrality `< 1e-8` stay unchanged.
- 2026-10-30 -- The `1e-8` target had no physical basis; it was a numerical
  convention inherited from the M1 tests. Physically meaningful targets are
  `μ` ~1e-3 eV (S15 uses 0.01%), density relative L2 ~1e-3, and pressure
  ~1%. The M1.5 spike pass criteria are therefore: no divergence, a
  perturbation response of order the perturbation (not amplified), and
  N-stability at the 1e-3 level, not `1e-8`.
- 2026-10-30 -- Decision F1: `average_atom` defaults to `N = 50`; the
  self-consistent N sweep uses `N = 30, 40, 50` with `max |Δμ| < 1e-2` eV.
  The fixed-potential `N = 50, 100, 150` convergence stays in the M1 test.
  High-precision eigensolves (B2) are deprioritized; the paper-aligned
  fallback is an M1.5 direct-G or resolvent reformulation.
- 2026-10-29 -- Initial `V_eff` is `-Z/r + Z/R` inside the cell, zero outside,
  so `V_eff(R) = 0` from the start.
- 2026-10-29 -- Free-electron subtraction is validated at fixed potential
  first, then in the self-consistent loop. It is a separate M2 sub-increment.
- 2026-10-29 -- New dependencies are allowed when they materially simplify a
  robust implementation, after a brief source/API check. Well-known/de facto
  packages such as `QuadGK` may be added without separate confirmation. Use
  `Pkg` to update `Project.toml`/`Manifest.toml`; do not hand-edit them.

## Notes

- S15 explicitly names KSDT (Karasiev et al., PRL 112, 076403 (2014)).
  SS23's text only writes `V_xc` and finite-temperature DFT, so the KSDT
  attribution belongs to S15 and the Tartarus lineage, not a direct SS23
  equation label.
- The M1 contour + Matsubara density remains the production density method.
  Free-electron subtraction changes the angular-momentum representation, not
  the contour formula.

## Blocker

2026-10-30 — M2a does not reach the required `1e-8` self-consistency. The M1
density layer is not a smooth function of `V_eff` at that scale, so no mixing,
iteration count, or acceleration can drive `max |ΔV_eff|` below ~1e-5. This is
the spec stop condition "self-consistent loop does not converge for the test
case with mixing in 0.3-0.5".

### What was tried

- Prototype SCF under `scratch/` using the M1 `fixed_potential_density` layer,
  a spectral cumulative Hartree integral (validated against an analytic density
  to ~3e-13), and a DVR cardinal interpolation of the smooth part of `V_eff`.
- Case `Z=5, A=10, rho_m=10, T=0.2` eV (`R=1.3882` Bohr), `N=40`, `l_max=4`,
  mixing 0.1/0.2/0.3/0.4, `max_iterations=300`: the residual stalls at
  3e-6 ... 1e-5 and oscillates. It never reaches `1e-8`.
- `l_max = 0, 1, 2` alone still floors at 2.5e-7, 5.8e-6, 1.4e-5.
- Milder cases (`Z=1,2`, `A=1,4`, `rho_m=10`, `T=0.2` eV) floor at
  4e-5 ... 4e-4.
- Larger `N` hits `ArgumentError("could not bracket the chemical potential")`
  in the M1 solver (`N=60`, `l_max=0`).

### Cause

The contour + Matsubara density sums over all `2N+l` Siegert pseudostates.
The high-|k| pseudostates have normalization constants `|N_n|` as small as
1e-11 (near-null SPS norms), so their normalized wavefunctions are huge
(`|P|^2` up to ~1e11 at `r ~ R`) and cancel catastrophically. The resulting
density carries a systematic Float64 error of ~1e-5 relative to `max |n_e|`.
That error is invisible in the M1 contour-invariance tests because the states
are fixed there, but it changes when `V_eff` changes.

Evidence at the `N=40`, `l_max=4` near-fixed point:

- Perturbing `V_eff` by `1e-8` changes `max |n_e|` by 4.8e-4 (relative 7.7e-6)
  and `μ` by 3e-5, independent of the perturbation size down to `1e-10`.
- Restricting the sum to `|k| <= 3` removes the noise (1.4e-8), but `|k| <= 3`
  changes the density value by 4% and `|k| <= 5` by 33%, so the high-k states
  are physically required, not negligible.
- The bound-state part is smooth (6e-9 for a `1e-8` perturbation); the contour
  and Matsubara parts carry the noise.
- Tighter `mu_tolerance` (`1e-12`, `1e-14`) does not help.
- Freezing the contour (`E_max`, `E_I`) does not help.
- LAPACK balancing before `geev` reduces the noise only from 4.8e-4 to 2.8e-4.
- Refining the contour quadrature (`n_horiz = n_return = 80`) does not help.

### Options (one decision needed before continuing)

A. Relax the M2 acceptance criteria to what the M1 layer can deliver, e.g.
   `max |ΔV_eff| < 1e-5` Ha and `|V_eff(R)| < 1e-5` Ha, while keeping the
   neutrality tolerance at `1e-8` relative (which is already met, ~1e-12).
   This changes the acceptance criteria in this SPEC.

### Resolution

2026-10-30 -- user chose option A. The self-consistency threshold becomes
`1e-4` Ha; the `V_eff(R) = 0` and neutrality thresholds stay at `1e-8`.
Work continues with M2a/M2b under the relaxed self-consistency threshold.
### 2026-10-30 update -- N-dependence also blocks criterion 5

After implementing M2a/M2b, the M1 density noise was measured against `N`. It
grows with `N`, so the self-consistent loop is only feasible for small bases:

- `Z=5, A=10, rho_m=10, T=0.2` eV, `l_max=4`, `mixing=0.4`, 100 iterations:
  - `N=40,50`: converge (`residV` ~ 4e-5 ... 1e-4).
  - `N=60`: converges with `l_max=4` (41 iterations) but not with
    `l_max=6,10`.
  - `N=70,80`: do not converge; `residV` stays ~5e-4 ... 3e-2.
  - `N=100,150`: wander and often throw `could not bracket the chemical
    potential`. `N=100, l_max=0` even ran to `μ = 362` Ha.
- `Z=1, A=4, rho_m=10, T=0.2` eV, `l_max=4`: `N=50` converges; `N=100`
  throws; `N=150` does not converge.
- `N=50, l_max=10` (the SPEC default `l_max`) converges in 20 iterations, so
  the spec's example case is feasible at `N=50`.
- The default `average_atom` value `N=150` therefore does not converge.

Self-consistent `N = 30, 40, 50` (Z=5, `l_max=10`) gives:

- density relative L2 changes: 2.8e-5 and 1.6e-4 (threshold 1e-2: pass),
- `max |Δμ| = 6.7e-3` eV (threshold `1e-4` eV: fails).

So criterion 5 as written (self-consistent `N=50,100,150`, `max|Δμ| < 1e-4`
eV) cannot be met by the current M1 layer, regardless of the self-consistency
tolerance change. The `μ` threshold is ~67x below the observed M1 noise floor
(~1e-4 Ha = 2.7e-3 eV).

Options (need one decision):

- F1 (recommended): Change the `average_atom` default `N` to 50. Make
  criterion 5 a self-consistent sweep over feasible `N = 30, 40, 50` and
  relax `max|Δμ|` to `1e-2` eV (density L2 threshold stays `1e-2`). Report
  the actual numbers in the test output.
- F2: Keep criterion 5's `N=50,100,150` but verify it on the fixed-potential
  M1 layer (already passing), and rely on criterion 1 for the self-consistent
  test. Change the `average_atom` default `N` to 50.
- F3: Fix the M1 density precision (option B) or implement M2c first
  (option C), then keep criterion 5 as originally written
  (`N=50,100,150`, self-consistent). More work, and M2c may not address the
  high-`k` conditioning.

### Resolution and revised spike criteria

2026-10-30 -- user chose F1. The `average_atom` default becomes `N = 50`, and
criterion 5 uses `N = 30, 40, 50` with `max|Δμ| < 1e-2` eV. The user also
approved replacing the unfounded `1e-8` spike target with physical criteria:

- no divergence or wandering to a spurious solution,
- density response to a `δV` perturbation is `O(δ)` (not amplified by the
  near-null states),
- `N = 50, 100, 150` agree at the physically meaningful level (`μ` ~1e-3 eV,
  density relative L2 ~1e-3).

B2 (high-precision eigensolve) is deprioritized. The next M1-level
investigation, if needed, is the S15-style direct Green function or a Siegert
resolvent formulation.

B. Fix the M1 density precision first (outside the current M2 scope): evaluate
   the contour/Matsubara sums with higher-precision accumulation (BigFloat or
   double-double) to recover the cancelled digits, or reformulate the density
   through the resolvent instead of the eigen-sum.

C. Implement M2c free-electron subtraction first and re-measure. The
   subtracted Green function may be much better conditioned; this is not
   guaranteed and would reorder M2a/M2c.

D. Change the density representation more broadly (not recommended silently).

### Current tree state

No tracked source has been modified. `scratch/` contains only the untracked
diagnostic scripts used above; it can be deleted or kept by the user.

## Implementation status

2026-10-30 -- M2a/M2b/M2c are implemented in `src/average_atom.jl` and tested
in `test/test_average_atom.jl`. Evidence from `nix develop`:

- Criterion 1: `average_atom(5, 10, 10, 0.2; N=50, l_max=10)` gives
  `converged = true`, 20 iterations, `self_consistency_residual = 9.89e-5`,
  and `|neutrality_residual|/Z = 3.4e-12`.
- Criterion 2: `V_eff(R) = -3.05e-11` Ha and `minimum(density) = 0.27`.
- Criterion 3: a constant `V_xc` leaves the result unchanged to `< 1e-3`; a
  linear `V_xc` changes `V_eff` (max ~0.11 Ha), keeps `V_eff(R) < 1e-8`, and
  still converges.
- Criterion 4: the unit round trip `average_atom` vs `_average_atom_bohr`
  matches to machine precision; the constants are tested directly.
- Criterion 5: self-consistent `N = 30, 40, 50` gives
  `max |Δμ| = 6.67e-3` eV and density relative L2 changes 2.8e-5 and 1.6e-4
  (thresholds 1e-2 and 1e-2).
- Criterion 6: at fixed `V = -Z/r + Z/R`, `N = 40`, subtraction at
  `l_max = 10` vs no-subtraction at `l_max = 20` gives relative L2 1.16e-7;
  on/off at `l_max = 20` gives 1.16e-7 (threshold 1e-6).
- Criterion 7: `just test` passes.
- Criterion 8: `treefmt --no-cache --fail-on-change src test` reports
  0 changed files.

Known limitation: `free_electron_subtraction = true` in the self-consistent
loop converges for `l_max = 4, 6` but not at `l_max = 10` for `Z = 5`
(`residV ~ 1e-3` after 100 iterations). The fixed-potential criterion passes;
the flag is off by default.

## Confirmed

2026-10-29 -- user: "yes"

2026-10-30 -- user approved the physical tolerances (1e-3 level, no 1e-8)
and asked to commit the M2 implementation with the documented limitations.
