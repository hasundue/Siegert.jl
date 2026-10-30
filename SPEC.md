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

- `max_r |V_eff,next(r) - V_eff(r)| < 1e-8` Ha, and
- `|4π ∫_0^R r^2 n_e dr - Z| / Z < 1e-8`.

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

- `average_atom(Z, A, rho_m, T; N=150, l_max=10, xc=nothing, mixing=0.4,
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
      `self_consistency_residual < 1e-8`, and
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

- [ ] `N = 50, 100, 150` on a self-consistent test case gives stable `μ` and
      density. Report the actual numbers in the test output. Initial thresholds:
      `max |Δμ| < 1e-4` eV and density relative `L2` change `< 1e-2`. If these
      are too tight for a legitimate case, stop and report rather than loosen
      silently.
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

## Confirmed

2026-10-29 -- user: "yes"
