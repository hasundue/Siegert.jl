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

## Confirmed

(awaiting explicit user yes)
