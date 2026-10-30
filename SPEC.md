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

Remaining B variant: biorthogonal/left-eigenvector normalization. The idea is
to pair each right eigenvector with its left partner and use the
biorthogonal denominator (which may be `O(1)`) instead of the tiny self-norm.
This is research-level (1-2 weeks), because the physical Green's function must
be re-derived in the extended space and validated against the existing M1
tests; the physical projection is still tiny, so success is not guaranteed.

Options:

- A. Accept the `N ~ 50` limitation and close M2.5 with this diagnosis (M2 as
  committed, plus the contour-panel note).
- B. B2 (biorthogonal/left-eigenvector normalization) as a bounded research
  spike. Expected 1-2 weeks; uncertain.
- C. M1.5: S15-style direct regular/irregular Green function or a Siegert
  resolvent, which avoids the near-null eigenvector normalization entirely.
- D. Proceed with M2.5b/M2.5c anyway (expected not to fix the root cause).

Recommendation: B2 only if the user wants to keep the Siegert spectral sum;
otherwise C is the more reliable reformulation. A remains valid given the
`N ~ 50` practical range.

## Confirmed

(awaiting explicit user yes)
