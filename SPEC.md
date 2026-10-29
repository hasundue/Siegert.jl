# SPEC -- Siegert Green's function and fixed-potential AA density (M1)

## Goal

Add the Siegert-state normalization, Green's function, DOS, and contour +
Matsubara density layer needed by the average atom model of Starrett & Shaffer,
*Phys. Rev. E* **107**, 025204 (2023). Given a fixed finite-range effective
potential `V_eff(r)` and ion-sphere radius `R`, compute `n_e(r)` and the
chemical potential `μ` that enforces charge neutrality. The continuum density
is computed only with the paper's contour + Matsubara formula; no second
continuum integrator is added.

This is milestone M1. The self-consistent AA loop (M2) and KSDT exchange-
correlation / pressure (M3) are out of scope.

## References

- C. E. Starrett and N. R. Shaffer, "Average Atom Model with Siegert States",
  *Phys. Rev. E* **107**, 025204 (2023), arXiv:2212.06229.
  Cited below as SS23.
- C. E. Starrett, "A Green's function quantum average atom model",
  *High Energy Density Physics* **16**, 18-22 (2015),
  doi:10.1016/j.hedp.2015.05.001. Cited below as S15. The contour and
  Matsubara construction is in S15 Sec. III and Appendix B.
- D. A. Batishchev and O. I. Tolstikhin, "Siegert pseudostate formulation of
  scattering theory: Nonzero angular momenta in the one-channel case",
  *Phys. Rev. A* **75**, 062704 (2007). Eq. (37) is the normalization used in
  SS23 Eq. (6); the reverse-Bessel zeros are already implemented.

## Current state

- `src/dvr.jl` builds the Gauss-Radau-Jacobi DVR basis with `α = 0`, `β = 2l`.
  The DVR basis matrix is diagonal: `π_i(x_j) = δ_ij / κ_i`.
- `src/sps.jl` assembles the TON98/BT07 matrices from the DVR basis and solves
  the linearized QEP. `SPS` returns `momenta`, raw `coefficients`, `grid`
  (radii), `hamiltonian`, `metric`, and `boundary`.
- `src/spectra.jl` implements the product formula for `S_l(k)` and the phase
  shift.
- Existing tests: `test/test_dvr.jl`, `test/test_sps.jl`.

A numerical check on the square well used by the existing tests shows that the
raw eigenvectors obey the SS23 Eq. (6) orthogonality with the DVR metric `ρ`
plus the surface term, within finite-basis accuracy, but they are not scaled to
the right-hand side `δ_nm`. Normalization is therefore required before any
Green's-function use.

## Scope (M1)

### 1. Normalize the Siegert states

For every state `n` of an `SPS`, compute

```text
N_n = (a/2) c_n^T ρ c_n
      + i P_n(R)^2 / (2 k_n)
        * [1 + Σ_{p=1}^l z_lp / (i k_n R + z_lp)^2]
```

where `c_n` is the raw DVR coefficient vector, `ρ` is the DVR metric,
`P_n(R) = a π_N^T c_n` is the boundary value of the radial function
`P_{n,l}(r) = r Σ_i c_{i,n} π_i(x(r))`, and `z_lp` are the reverse-Bessel
zeros. The factors `a/2` and `a` come from `dr = (a/2) dx` and
`P = r φ`; the earlier draft omitted them. This is SS23 Eq. (6) with `m = n`
and BT07 Eq. (37). Rescale `c_n ← c_n / sqrt(N_n)`. The square root uses
Julia's principal branch; the overall sign cancels in every M1 observable
because they are quadratic in `P_n`.

Store the raw normalization constants and the DVR object in `SPS`, and make
`SPS.coefficients` the normalized coefficients. Existing SPS tests that use
only `momenta` must continue to pass.

### 2. Evaluate wavefunctions on a common radial grid

The DVR grids depend on `l`, so states with different `l` cannot be summed on
their native grids. Use the `l = 0` DVR grid as the common radial grid for all
partial waves. Add an evaluator for the DVR/FBR basis at arbitrary `x`:

```text
φ_n(x) = sqrt(w(x) / h_n) P_n^{(0, 2l)}(x)
π_i(x) = Σ_n T_{ni} φ_n(x)
P_{n,l}(r) = r Σ_i c_{i,n} π_i(x(r))
```

using the existing `T`, `w`, `h_n`, and `jacobi` machinery. Evaluate every
`l` state on the common grid; precompute the matrix `P[grid, state]` once per
`l` and reuse it for all contour and Matsubara points.

The `l = 0` Gauss-Radau weights provide the radial quadrature. For
`r = R(1+x)/2`, the volume integral is

```text
∫_0^R r^2 f(r) dr = (R^3/8) Σ_i ω_i (1+x_i)^2 f(r_i).
```

### 3. Green's function and DOS

Implement the partial-wave Green's function diagonal on the common grid

```text
G_l(r,r;k) = Σ_n P_{n,l}(r)^2 / [k_n (k - k_n)].
```

`k` may be complex; the sum runs over all `2N+l` Siegert states. Do not
complex-conjugate `P_n`.

Implement the partial and total DOS as

```text
χ_l(k) = -(1/π) Im ∫_0^R dr 2(2l+1) G_l(r,r;k)
χ(k)   = Σ_l χ_l(k).
```

This is the published SS23 Eqs. (16)-(17). The DOS is a volume integral:
the 3D volume element `4π r^2 dr` cancels the `1/r^2` in the local density
`G_l/r^2`, so no `/r^2` appears. Equivalently, the 3D Green's function is
`G_3D = Σ_l (2l+1)/(4π) G_l/r^2`, and the DOS is
`-(2/π) Im ∫_0^R 4π r^2 dr G_3D`. The single-state formulas SS23
Eqs. (18)-(19) follow from the same cancellation and the normalization
`∫_0^R P_n^2 dr = 1`. The earlier SPEC draft incorrectly treated the
missing `/r^2` as a typo; the published PRE PDF resolves this.

### 4. Electron density

Use the decomposition

```text
n_e(r) = n_bound(r) + n_cont(r).
```

Bound states are the Siegert states with `Re(k_n) = 0` (within numerical
tolerance) and `Im(k_n) > 0`. Their contribution is

```text
n_bound(r) = Σ_{l} Σ_{n bound}
    f(ε_n, μ) 2(2l+1)/(4π) P_{n,l}(r)^2 / r^2,
ε_n = k_n^2 / 2.
```

The continuum part uses SS23 Eq. (13) with `E_min = 0` by default:

```text
n_cont(r) = (1/π) Im ∫_C dε f(ε, μ)
                Σ_l 2(2l+1)/(4π) G_l(r,r;k) / r^2
          + 2 T Re Σ_{j=1}^{N_mat}
                Σ_l 2(2l+1)/(4π) G_l(r,r;k_j) / r^2,
```

with `k = sqrt(2ε)`, `k_j = sqrt(2 z_j)`, and Matsubara energies
`z_j = μ + iπ(2j-1)T`. `E_min = 0` is the simplest choice from S15 Appendix B:
all bound states are added explicitly, so there is no double counting. Support
a negative `E_min` as an option, in which case only bound states below `E_min`
are added explicitly.

Contour construction (S15 Appendix B, Fig. 1): the rectangular path

```text
E_max → E_max + i E_I → E_min + i E_I → E_min
```

traversed from `E_max` to `E_min`, with the real-axis segment closing the
rectangle. Use 4-point Gauss-Legendre per panel: 1 panel on the short
`E_max → E_max + i E_I` section, 20 on the horizontal section, and 20 on the
`E_min + i E_I → E_min` section (164 points total by default). The `E_max`
choice is the energy at which the Fermi factor has decayed to `1e-10`,
computed from the current `μ`; clamp it to at least 0. Choose `E_I = 1` Ha
when there are no Matsubara poles inside the rectangle. Otherwise choose
`N_mat` so that `E_I = 2 N_mat π T` is as close to 1 Ha as possible while
remaining above the highest enclosed pole; count only poles with
`E_min < Re(z_j) < E_max`. This is the S15 scheme.

Implement the contour integral literally as written. If a sign or prefactor is
wrong, the bound-state limit, free-particle limit, or neutrality test will
fail; stop and report instead of fitting a factor.

### 5. Chemical potential

Solve

```text
Z = 4π ∫_0^R r^2 n_e(r) dr
```

for `μ` by bisection with an expanding bracket. Recompute the contour
parameters `E_max`, `E_I`, and `N_mat` for each trial `μ`. Require `T > 0`;
for `T → 0`, use a small positive temperature in tests rather than a special
zero-temperature path. Return the neutrality residual and the final contour
parameters in the result.

### 6. Fixed-potential driver

Implement a driver that

1. builds `SPS(N, l, V, R)` for `l = 0:l_max`;
2. chooses the `l = 0` grid and weights as the common quadrature;
3. precomputes all normalized wavefunctions on the common grid;
4. solves for `μ` and `n_e`;
5. returns the density, `μ`, states, grid, contour metadata, and residual.

This driver does not update `V`. It is the M1 deliverable and the input to the
M2 self-consistent loop.

## API

Proposed public names (one export per line, as required by the project style):

- `SPS`: gains `normalizations::Vector{ComplexF64}` and `dvr::DVR`; its
  `coefficients` become normalized.
- `siegert_wavefunction(sps::SPS)` -> `Matrix{ComplexF64}` on `sps.grid`.
- `siegert_wavefunction(sps::SPS, r::AbstractVector)` -> `Matrix{ComplexF64}`
  on the requested radii.
- `siegert_green(sps::SPS, k::Number; r=sps.grid)` -> `Vector{ComplexF64}` of
  diagonal Green's-function values on `r`.
- `siegert_dos(sps::SPS, k::Number; r, weights)` -> `Float64` for one `l`.
- `siegert_dos(states::AbstractVector{SPS}, k::Number; r, weights)` ->
  `Float64` for the total.
- `fermi_dirac(ε, μ, T)` -> `Float64`.
- `matsubara_energies(μ, T; E_min, E_max, E_I)` -> `Vector{ComplexF64}`.
- `siegert_density(states, μ, T; ...)` -> `Vector{Float64}` on the common grid.
- `chemical_potential(states, Z, R, T; ...)` -> `(μ, density, residual)`.
- `fixed_potential_density(V, Z, R, T; N=150, l_max=10, ...)` ->
  `AADensity`.
- `AADensity`: result struct with `density`, `chemical_potential`, `states`,
  `grid`, `radius`, `temperature`, `neutrality_residual`, and `contour`.
- `Contour`, `contour_parameters`, `common_radial_grid`, and `volume_weights`
  are exported helpers so the driver result and the quadrature are testable
  without private functions.

Names may be adjusted during implementation, but the driver and result must
remain testable without reaching into private functions.

## Units and conventions

- Internal calculations use Hartree atomic units: `ε = k^2/2`, `T` in Hartree,
  `R` in Bohr.
- M1 takes `T` in Hartree. The M2 public `average_atom(Z, A, ρ_m, T)` will
  accept `T` in eV and `ρ_m` in g/cm³ and convert at the boundary.
- `V_eff(r)` is a function zero for `r ≥ R`; the DVR grid avoids `r = 0`.
- The density includes spin degeneracy 2 and the `(2l+1)` orbital degeneracy
  through the `2(2l+1)/(4π)` factor.

## Acceptance criteria

- [ ] `SPS` states satisfy SS23 Eq. (6) for all states of a square well with
      `N = 12`, `l = 0, 1, 2`: the normalized inner product is `δ_nm` to
      `1e-8` or better, including complex states.
      verify: `just test average_atom`

- [ ] For bound states of the existing square well, the ordinary normalization
      `∫_0^∞ P_n(r)^2 dr = 1` holds to `1e-8` or better when the DVR volume
      integral is combined with the analytic tail `P_n(R)^2/(2κ)`.
      verify: `just test average_atom`

- [ ] `siegert_wavefunction(sps, r)` reproduces the native-grid values on
      `sps.grid` to `1e-10` for `l = 0, 1, 2`, and evaluates all `l` on the
      `l = 0` grid without errors.
      verify: `just test average_atom`

- [ ] For a bound state, the residue of `siegert_green` at `k → k_n` matches
      `P_n(r)^2 / k_n` on the common grid to a relative tolerance of `1e-6`.
      verify: `just test average_atom`

- [ ] For the free potential `V = 0`, the DOS `siegert_dos` is positive for
      real `k > 0` and approaches the free-electron value `V k / π^2` as `N`
      and `R` grow. The test uses a tolerance and a convergence check rather
      than a hard absolute value.
      verify: `just test average_atom`

- [ ] With a deep square well and `μ` chosen so the continuum is empty, the
      contour + Matsubara continuum density is negligible and the total
      density equals the bound-state sum to `1e-10` relative tolerance.
      verify: `just test average_atom`

- [ ] On a smooth finite-range test potential with continuum occupation, the
      density is invariant under the contour choices: doubling the panel
      counts, changing `E_I` within the allowed range, and using `E_min = 0`
      versus a negative `E_min` (with the corresponding bound states moved to
      the explicit sum) agree to `1e-6` relative tolerance.
      verify: `just test average_atom`

- [ ] `fixed_potential_density` enforces neutrality:
      `|4π ∫_0^R r^2 n_e dr - Z| / Z < 1e-8`, the density is non-negative
      within roundoff, and `μ` is finite for a test case with both bound and
      continuum electrons.
      verify: `just test average_atom`

- [ ] `N = 50, 100, 150` on the same test potential gives a stable `μ`
      (change below `1e-4` Ha) and density (relative `L2` change below
      `1e-2`); report the actual numbers in the test output.
      verify: `just test average_atom`

- [ ] Existing DVR and SPS tests still pass.
      verify: `just test dvr sps` or `just test`

- [ ] Formatting is clean.
      verify: `treefmt -c`

## Out of scope

- The self-consistent AA loop: updating `V_eff` from `n_e`, the Hartree term,
  and `V_xc[n] - V_xc[n(R)]`.
- All exchange-correlation models, including finite-temperature LDA/KSDT.
- Mass density / atomic mass to ion-sphere radius conversion, and the public
  `average_atom(Z, A, ρ_m, T)` API.
- Pressure, entropy, free energy, equation of state, and comparison with
  Tartarus.
- Free-electron subtraction / acceleration of the angular-momentum sum.
- Off-diagonal Green's function `G_l(r,r';k)` for `r ≠ r'`.
- Relativistic effects, non-integer `l`, multi-channel scattering, and
  non-spherical potentials.
- Reproduction of the SS23 figures for Al, Fe, C, or Cu.
- Fixing unrelated stale examples or docstrings (for example
  `examples/basis.jl` or the `square_well` docstring).
- New dependencies. The 4-point Gauss-Legendre panels are implemented
  directly; do not add `QuadGK` or another quadrature package for M1.

## Stop conditions

- Eq. (6) cannot be made to hold for all states of the `N = 12` square-well
  test after two normalization attempts -> write the cause under `## Blocker`,
  stop.
- The contour + Matsubara density disagrees with the bound-state limit or with
  contour invariance after two attempts -> write the cause, stop. Do not add a
  second continuum integrator to hide the disagreement.
- The sign or prefactor of SS23 Eq. (13) cannot be resolved from the
  free-particle and neutrality tests -> stop and report.
- Neutrality fails because bound-state wavefunction tails outside `R` are not
  counted in the cell integral -> stop and report. Do not add an ad hoc tail
  correction without a decision.
- The DOS or density normalization becomes material: for example, a factor of
  `4π` or `r^2` cannot be resolved from the volume-integral derivation and the
  published equations -> stop and report. Do not silently change the formula
  or fit a factor.
- A required reference cannot be obtained or opened -> stop and report which
  reference is missing. Do not substitute a secondary summary or an
  unauthorized copy; the user will provide the published PDF.
- A change is needed outside `src/dvr.jl`, `src/sps.jl`, a new
  `src/average_atom.jl` (or `src/greens.jl`), `src/Siegert.jl`,
  `test/test_average_atom.jl`, and the related docstrings -> stop.
- A new dependency is needed -> stop before editing `Project.toml` or
  `Manifest.toml`.
- Existing SPS tests fail -> stop and report the discrepancy.

## Decisions

- 2026-10-29 -- M1 is the Siegert normalization / Green's function /
  fixed-potential density layer; self-consistency, XC, and pressure are later.
- 2026-10-29 -- Contour + Matsubara is the only production continuum density
  method from M1 onward. Validation uses analytic limits, bound-state limits,
  contour invariance, and neutrality; it does not introduce a second
  continuum integrator.
- 2026-10-29 -- `E_min = 0` by default, with all bound states added explicitly.
  This is the simplest choice in S15 Appendix B and avoids double counting.
  Negative `E_min` is an option, not the default.
- 2026-10-29 -- Use 4-point Gauss-Legendre panels as in S15 Appendix B; no new
  quadrature dependency.
- 2026-10-29 -- DOS uses the published volume-integral form:
  `χ_l(k) = -(1/π) Im ∫_0^R dr 2(2l+1) G_l(r,r;k)`, with no `/r^2`. The
  `4π r^2` volume element cancels the `1/r^2` in `G_3D`; the single-state
  Eqs. (18)-(19) are consistent with `∫P_n^2 dr = 1`. The published PRE PDF
  (Zotero attachment `YZ7IGZKR`) and the volume-element derivation agree.
- 2026-10-29 -- Treat normalization and prefactors as high-risk: TON98, BT07,
  and SS23 are not fully consistent in some related formulas. If a factor
  discrepancy becomes material, stop and ask for the published PDF/erratum
  instead of fitting the result.
- 2026-10-29 -- If a needed reference cannot be accessed, stop and report it;
  the user will provide the published PDF or an erratum. Do not substitute a
  secondary summary or an unauthorized copy.
- 2026-10-29 -- The cell density and DOS use the paper's domain `0 ≤ r ≤ R`.
  Bound-state tails outside `R` are not counted. If this prevents neutrality,
  stop and report; do not patch the density with an unagreed tail term.
- 2026-10-29 -- Normalize inside `SPS`; `SPS.coefficients` become normalized,
  and `SPS` stores `normalizations` and `dvr`.
- 2026-10-29 -- M1 uses Hartree units. M2 will convert eV and g/cm³ at the
  public API boundary.
- 2026-10-29 -- M2/M3 XC remains pluggable: a simple XC first, KSDT later.
  This is a roadmap decision, not part of M1.
- 2026-10-29 -- Implementation: the radial Siegert function is
  `P = r φ`, so the volume term is `(a/2) c^T ρ c` and the boundary value is
  `P(R) = a π_N^T c`. With these factors SS23 Eq. (6) holds to about `1e-12`
  for the square-well test in `l = 0, 1, 2`.
- 2026-10-29 -- Implementation: the Jacobi norms use `loggamma` and the
  Jacobi polynomials use a local three-term recurrence. The external
  `jacobi`/`gamma` path overflows near `N = 100`, which would block the
  required `N = 150` convergence test. Existing DVR and SPS tests are
  unchanged.
- 2026-10-29 -- Implementation: the free-electron DOS test uses `N = 5, 6, 7`
  and `R = 1, 2`. For `V = 0`, some SPS bilinear norms approach zero at
  `N >= 8`, so the real-axis DOS becomes ill-conditioned; the contour +
  Matsubara density is unaffected because it stays away from the real axis.

## Roadmap (not accepted by this SPEC)

- **M2**: self-consistent AA loop. Build `V_eff` from `-Z/r`, Hartree, and a
  pluggable `V_xc`; iterate to convergence; add free-electron subtraction for
  the angular-momentum sum; expose `average_atom(Z, A, ρ_m, T)`.
- **M3**: finite-temperature KSDT exchange-correlation, pressure/entropy,
  and comparison with the SS23 figures and Tartarus.

## Resolved blocker

2026-10-29 -- The published PRE PDF is now in Zotero (attachment `YZ7IGZKR`,
DOI 10.1103/PhysRevE.107.025204, PDF metadata title "Average-atom model with
Siegert states"). It confirms that SS23 Eqs. (16)-(19) have no `/r^2`. This is
not a typo: the DOS is a volume integral and the `4π r^2` element cancels the
`1/r^2` in `G_3D`. The earlier `/r^2` concern was an error in the SPEC draft,
not in the paper. M1 is unblocked. No code has been changed.

## Confirmed

2026-10-29 -- user: "M1 を開始してください"
