# SPEC -- M2.5: stable density evaluation for larger Siegert bases

## Goal

Make the average-atom density and self-consistent loop stable up to `N = 150`
at physically meaningful tolerances, without changing the M1 public API
semantics and without fitting ad hoc corrections.

The work is ordered as agreed: **M2.5a contour/pole audit, M2.5b faithful
free-electron correction, M2.5c SCF mixing/acceleration**. Each sub-increment
stops and reports before the next one starts if it changes scope.

## Background

- M2 implementation: `src/average_atom.jl`, `src/sps.jl`, `src/dvr.jl`,
  `test/test_average_atom.jl`, `README.md`. M2 is committed with documented
  limitations (`N ~ 50`, free subtraction only fixed-potential-validated).
- SS23 uses TON97/TON98/BT07 and the same Jacobi-DVR basis and normalization;
  a different basis is ruled out as the explanation.
- SS23 (arXiv:2212.06229) uses the density-level free-electron acceleration of
  Tartarus (arXiv:1804.01613) Sec. 3.3; the published text only says
  `l > l_max` is treated as free electrons.
- Tartarus Sec. 3.3 uses `n_e = orbitals + ctm(l_max+1..l_con)
  - free(0..l_con) + n_e^0`, with the high-`l` tail computed by orbitals on a
  fixed `sqrt(epsilon)` grid, `l_con` from `TOL = 1e-4`, and `n_e^0` the
  uniform free-electron density. Tartarus also uses `alpha = 0.1` simple
  mixing and Eyert/Anderson acceleration.
- Physical tolerances from M2: self-consistency `< 1e-4` Ha, neutrality
  `< 1e-8` relative, `V_eff(R) = 0` within `1e-8`, `mu` N-stability `~1e-3`
  eV, density relative L2 `~1e-3`. No `1e-8` target.

## M2.5a: contour and Siegert-pole audit

Diagnose before changing anything.

1. For a converged M2 potential (e.g. `Z=5, A=10, rho_m=10, T=0.2` eV at
   `N=50`), list the Siegert eigenvalues inside the contour
   `0 < Re(epsilon) < E_max`, `0 < Im(epsilon) < E_I`, with their normalization
   constants and residues.
2. Check whether SS23 Eq. (13), which only includes Matsubara residues, is valid
   for that contour, and compare `E_min`, `E_max`, `E_I` with S15 Appendix B.
3. Reproduce the `N = 100` bracket garbage and check whether it correlates with
   an enclosed pole, with `E_max` approaching the top of the pseudospectrum, or
   with near-null cancellation.
4. If the contour is at fault, propose a contour rule (no fitted correction)
   and test `N = 50/100/150`.

Acceptance criteria:

- [ ] A written diagnosis with numbers for the enclosed poles and the
      `N = 100` garbage.
      verify: report in this SPEC and `just test average_atom`.
- [ ] If a contour change is made, `average_atom(5, 10, 10, 0.2; N=100/150,
      l_max=10)` converges at the physical tolerances.
      verify: `just test average_atom`
- [ ] If the contour is valid, record that and proceed to M2.5b without
      changing it.
      verify: report in this SPEC.

Stop conditions: the same failure twice; an enclosed pole cannot be removed
without a fitted factor; the diagnosis needs code outside the listed files.

## M2.5b: faithful free-electron correction

Implement the Tartarus/SS23 density-level correction:

`n_e = n_Siegert(l <= l_max) + n_ctm(l_max+1..l_con)
      - n_free(0..l_con) + n_e^0`

with the free reference analytic (Bessel/Hankel, S15 Eq. (14)) and the high-`l`
tail either free or from a stable orbital representation. Determine `l_con` by
the Tartarus `TOL = 1e-4` criterion. Validate at fixed potential first, then in
the self-consistent loop.

Acceptance criteria:

- [ ] Fixed-potential correction matches the no-correction density at a larger
      `l_max` to `1e-6` relative L2, and on/off agree at the same large
      `l_max`.
      verify: `just test average_atom`
- [ ] `average_atom(5, 10, 10, 0.2; N=100/150, l_max=10)` converges at the
      physical tolerances, or a documented reason is recorded.
      verify: `just test average_atom`

Stop conditions: the free reference cannot be defined consistently with the
finite cell; a new dependency would be needed; the correction requires a fitted
factor.

## M2.5c: SCF mixing/acceleration

Only if a/b leave a slow or oscillatory convergence.

- Test `alpha = 0.1` simple mixing first.
- Eyert/Anderson acceleration only after explicit agreement, because it changes
  the M2 SPEC decision (mixing 0.4, no accelerator).
- Also test a Thomas-Fermi initial guess only if agreed.

Acceptance criteria:

- [ ] `N = 100/150` converge within a stated iteration budget, with the mixing
      or accelerator documented.
      verify: `just test average_atom`

Stop conditions: a mixing schedule or accelerator is needed but not agreed.

## Out of scope

- M1.5 (S15 direct Green function or Siegert resolvent). Only reconsider if
  M2.5a-c fail.
- KSDT exchange-correlation, pressure, entropy, free energy, equation of state
  (M3).
- Changing the M1 public API semantics.
- Relativistic effects, non-spherical potentials, multi-channel scattering.

## Files

`src/average_atom.jl`, `src/sps.jl`, `src/dvr.jl`,
`test/test_average_atom.jl`, `README.md`, `SPEC.md`.

## Decisions

- 2026-10-30 -- M2 is closed with documented limitations (commit `7a8f16e`).
  M2.5 is ordered M2.5a -> M2.5b -> M2.5c.
- 2026-10-30 -- M1.5 (direct Green / resolvent) is deprioritized because the
  basis is the same as in SS23; pursue contour, free correction, and mixing
  first.
- 2026-10-30 -- The physical tolerances from M2 carry over; no `1e-8` target.

## M2.5a findings (2026-10-30)

The contour is not the primary cause. The high-`|k|` near-null Siegert states
are.

- At the converged `N=50, l_max=10` potential, the default contour encloses no
  upper-half Siegert pole. The nearest is `l=1, Re=0.680, Im=1.978`,
  `|norm|=2.0e-2`, outside the default `E_I=1.016`.
- Contour refinement: default 20 panels vs 80 panels changes the density by
  1.2e-5. With refined panels, `E_I=1.2` agrees with the base to 1.3e-8;
  `E_I=1.5` still differs by 5.1e-3 (the contour passes near the `Im=1.978`
  pole), and `E_I=2.0` by 7.6e-5. So the default contour is valid but has a
  ~1e-5 quadrature error; larger `E_I` is not a fix.
- The `N=100` failure at the wandered potential has `poles_in=0`, closest pole
  distance ~1.0, and contour refinement changes only 1e-4 ... 1e-3. Yet the
  density has `nmin ~ -1.3e3` and the neutrality residual ~-150. The contour
  is not producing that error.
- Evaluating the **same** converged `N=50` potential on larger grids (fixed
  physical potential, `mu=1.8466`):

  | N | nmin | min \|norm\| | near-null count (<1e-8) |
  | ---: | ---: | ---: | ---: |
  | 50 | +0.27 | 2.6e-13 | 629 |
  | 60 | -0.035 | 5.3e-14 | 827 |
  | 80 | -1.16 | 4.1e-15 | 1233 |
  | 100 | -80.1 | 5.3e-16 | 1649 |
  | 120 | -607 | 9.5e-17 | 2075 |

  The density degrades with `N` even when the potential is held fixed. This is
  a normalization/cancellation limit of the high-`|k|` Siegert states, not an
  SCF feedback effect.
- A `|k|` cutoff at `N=100` (same fixed potential): `|k| <= 3` is positive
  (`nmin=0.35`) but misses 6.4 electrons (`f=+6.4`); `|k| > 3` carries the
  negative garbage (`nmin=-82`). So truncation alone breaks neutrality, and the
  lost tail is not negligible.
- `free_electron_subtraction=true` leaves the `N=100` density unchanged
  (`nmin=-80.1`, `f=-1.88`): it subtracts the analytic free partials
  (`n_free ~ 0.1`), not the noisy high-`k` Siegert tail.

## Blocker / decision needed

M2.5b (the `l > l_max` free correction) and M2.5c (mixing/acceleration) cannot
fix a high-`k` normalization failure inside each `l <= l_max`. Proceeding with
them as written would not address the observed `N=100/150` garbage.

### B1 diagnostic (2026-10-30)

B1 was "use BT07 Eq. (36) in the DVR basis, or evaluate the normalization in
higher precision, to reduce the near-null amplification". A direct test rules
this out: for the smallest-norm state at `N=100` (same fixed `N=50` potential),
`l=10`, `k = 2.52e-9 - 645.46i`,

- `|norm_f64| = 5.263865e-16`, `|norm_big| = 5.263862e-16` from the same
  Float64 eigenvector,
- relative difference `6.1e-7`,
- `interior ~ 3.21e-6`, `|surface| ~ 1.0`, `|P| ~ 0.064`, and the two `3.2e-6`
  terms cancel to `5e-16` (about 10 digits).

The normalization *evaluation* is already accurate. The tiny norm is genuine,
and the loss comes from the Float64 eigenvector plus the intrinsic
cancellation. B1 (DVR form or higher-precision normalization only) is ruled
out.

### B2 diagnostic (2026-10-30)

The biorthogonal/left-eigenvector variant was tested on the `l=0` linearized
QEP at `N=100` with the fixed `N=50` potential.

- The standard left/right eigenvectors of the companion matrix are extremely
  ill-conditioned (`cond(Y) ~ 9.1e13`, `cond(Z) ~ 2.0e14`). The standard
  overlap `|z^T y|` is tiny (`~1.3e-13`), so a naive left/right pairing does
  not help.
- The correct QEP biorthogonal form uses the left QEP eigenvector `d_n` and
  `ov_n = d_n^T (B + 2 lambda_n I) c_n`. For the smallest-norm state
  (`lambda = 539.42`):
  - `|norm| = 8.26e-16`, `|ov| = 1.24e-8`, `|d^T c| = 3.44e-6`,
  - `|d^T c| / |ov| ~ 2.8e2`, versus the current `|c^T rho c| / |norm| ~ 8e9`.
- Pointwise residue at the grid points, `r^2 c_j d_j / (kappa_j^2 ov)` versus
  the current `r^2 c_j^2 / (kappa_j^2 norm)`:
  - current max ~`6.0e12`, median ~`1.2e-4`,
  - biorthogonal max ~`4.0e5`, median ~`5.1e-14`.

The biorthogonal pairing reduces the amplification by about `1.5e7` (from
`1e12` to `1e5`). This is a strong positive signal, but not a proof: the
physical Green's function must be re-derived with the metric, boundary, and
contour/Matsubara factors, and validated against the existing M1 tests. The
pointwise `4e5` still implies a few digits of cancellation, so success is
plausible but not guaranteed.

Remaining B variant: B2 (QEP biorthogonal/left-eigenvector normalization).

Options:

- A. Accept the `N ~ 50` limitation and close M2.5 with this diagnosis (M2 as
  committed, plus the contour-panel note).
- B. B2 as a bounded research spike (1-2 weeks). The diagnostic is promising;
  start with a fixed-potential `l=0` prototype at `N=100`, compare against the
  `N=50` result, and stop if the density is not positive and stable.
- C. M1.5: S15-style direct regular/irregular Green function or a Siegert
  resolvent, which avoids the near-null eigenvector normalization entirely.
- D. Proceed with M2.5b/M2.5c anyway (expected not to fix the root cause).

Recommendation: B2 is now the preferred next step for keeping the Siegert
spectral sum; C remains the fallback. A is still valid for closing M2.5.

### B2 prototype result (2026-10-30)

The first B2 prototype (replace the self-normalized `P_n^2 / N_n` by the QEP
biorthogonal `P_n^R P_n^L / ov_n`, with `ov_n = d_n^T (B + 2 lambda_n I) c_n`)
**failed**:

- At `N=50` with the fixed physical potential, the biorthogonal density has
  `nmax = 21.7` vs the current `92.1`, and `biorth/current` is not a constant
  (median ~0, mean ~5.4). The shape and scale are wrong, not just a factor.
- At `N=100` the biorthogonal density has `nmin = -30.4`, so it is not stable.

The reason is that the QEP biorthogonal residue is the residue of the
*coefficient-space QEP resolvent*, not of the physical Siegert Green's
function. The physical `G_l` uses the self-normalized SPS (the Siegert inner
product), and the metric `rho`/boundary terms are missing from the naive
`P^R P^L / ov` replacement. The promising `4e5` diagnostic was for the QEP
resolvent, so it does not transfer directly.

Remaining B2 variant: use the extended-space self-orthogonality
`y_n^T W y_m = 2 lambda_n delta_nm` with the correct weight `W` for the
project's companion matrix, and derive the physical projection with `rho` and
the boundary terms. This is additional research (uncertain), not a drop-in
normalization change.

Revised options:

- A. Accept the `N ~ 50` limitation and close M2.5 with this diagnosis.
- B2-W. Continue B2 with the correct `W`-weighted extended-space derivation.
  Research-level; success uncertain.
- C. M1.5: S15-style direct regular/irregular Green function or a Siegert
  resolvent, which avoids the near-null eigenvector normalization entirely.
- D. Proceed with M2.5b/M2.5c anyway (expected not to fix the root cause).

Recommendation: C is now the more reliable route. B2-W is worth continuing
only if keeping the spectral sum is a hard requirement. A remains valid.

### C prototype result (2026-10-30)

S15-style direct Green function implemented for `l=0` (RK4 radial
regular/irregular solutions, Juelich matching, S15 Eq. (14) free form):

- Free case (`V=0`): `direct_G` matches `r^2 j_l(kr) h^+_l(kr)` to `1e-9`.
- Interacting `N=50, l_max=0` potential: the direct continuum matches the
  current Siegert continuum with a *constant* prefactor `2.0000000007` over the
  grid. The factor `2` is the spin convention in the project's
  `2(2l+1)/(4π)`; with `pref=2` the two continua agree.
- Total density (Siegert bound part + direct continuum, `pref=2`):
  - `nsteps=2000`: `nmin=0.1466`, `nmax=92.1`, relative L2 vs current
    `1.5e-4`,
  - `nsteps=4000`: `nmin=0.1466`, `nmax=92.1`, relative L2 vs current
    `6.3e-5`, change vs `nsteps=2000` `9.1e-5`.

The direct method is independent of the Siegert basis size, so it avoids the
near-null problem. The naive implementation costs ~36 s per density at
`nsteps=4000` (heavy allocations), so the production backend needs
optimization (cached k-points, fewer allocations, coarser/adaptive grid,
possibly parallel contour points).

Next step for C: optimize the `l=0` density, extend to `l>0`, and integrate it
as the density backend while keeping the Siegert states for the bound sum, DOS,
and resonance API. Then run a short self-consistent check at `N=100/150` with
the physical tolerances.

### C factor and l>0 validation (2026-10-30)

- Free `l=0,1,2` checks: `direct_G` matches `r^2 j_l(kr) h^+_l(kr)` with a
  constant ratio `1` (spread across `r` ~1e-6 ... 1e-8). The direct Green
  function itself is l-independent; the factor `2` is not in `direct_G`.
- The prefactor `2` is the structural ratio between the SPS full-spectrum
  Green function (`sum_{n=1}^{2N+l} P_n^2 / (k_n(k-k_n))`, both `+k` and `-k*`
  branches) and the S15 outgoing/Juelich representation. It is not a fitted
  shape correction. The free `N=6` ratio is ~1.9 (finite-N), and the interacting
  `N=50` ratio is `2.0000000007`.
- The factor must be a documented module constant with tests (free density,
  bound limit, l-independence), not a per-case calibration.

### Optimization status (2026-10-30)

A preallocated RK4 with V precomputed on a uniform fine grid is ~400x faster
(0.05 s vs ~20 s per density), but it is not yet accurate: near `r=1e-6` the
Coulomb `-Z/r` varies too rapidly for the midpoint V interpolation, giving a
~5% error at a single `k` and relative L2 ~3.5 in the density. The validated
slow version uses the exact V at every RK4 stage. Next optimization: a
nonuniform grid clustered at the origin (or exact V at stage points), then
re-validate against the `N=50` density before integrating.

## Confirmed

(awaiting explicit user yes)
