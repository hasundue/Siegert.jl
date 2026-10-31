# SPEC -- M3: finite-temperature XC and equation-of-state observables

## Goal

Add a finite-temperature exchange-correlation (XC) model and thermodynamic
observables to the average-atom (AA) layer, using the existing direct-G density
backend as the default. Keep the M2 physical tolerances and do not fit ad hoc
corrections. Reproduce the physical content of SS23's DOS and pressure results
and make the model usable for low-temperature liquid-metal cases.

The work is ordered as agreed: **M3a corrKSDT XC -> M3b thermodynamic
observables/API -> M3c validation**. Each stage stops and reports before the
next starts if it changes scope.

## Background

- M2/M2.5 implementation: `src/average_atom.jl`, `src/sps.jl`, `src/dvr.jl`,
  `test/test_average_atom.jl`, `README.md`.
- `average_atom` returns `AverageAtomResult`; the default density backend is
  `:direct_green`, which is stable at `N = 100/150`. `:siegert` remains for the
  M1 spectral path, DOS, and resonances.
- The current `xc` is a callable `xc(n, r) -> V_xc`; the driver applies the
  `V_xc[n(R)]` subtraction. There is no `f_xc`, so the result has no XC free
  energy, entropy, or pressure.
- M2 tolerances: self-consistency `< 1e-4` Ha, neutrality `< 1e-8` relative,
  `V_eff(R) = 0` within `1e-8`, `mu_tolerance = 1e-10` Ha absolute. Do not chase
  a `1e-8` potential residual.
- The Siegert spectral density is unreliable above `N ~ 50` (near-null
  high-`|k|` states); SS23's `N = 150` N-convergence figures (Figs. 7/8) are a
  separate numerical-method limitation. M3 uses the direct-G density and does
  not reopen that problem.
- The physical target is SS23's DOS and pressure results, not the exact
  N-convergence of the Siegert spectral density.

### References

Core:

- [SS23] C. E. Starrett and N. R. Shaffer, "Average Atom Model with Siegert
  States," Phys. Rev. E **107**, 025204 (2023),
  doi:10.1103/PhysRevE.107.025204, arXiv:2212.06229.
- [S15] C. E. Starrett, "A Green's function quantum average atom model," High
  Energy Density Phys. **16**, 18-22 (2015),
  doi:10.1016/j.hedp.2015.05.001.
- [Tartarus] C. E. Starrett, N. M. Gill, T. Sjostrom, and C. W. Greeff, "Wide
  ranging equation of state with Tartarus: A hybrid Green's function/orbital
  based average atom code," Comput. Phys. Commun. **235**, 50-62 (2019),
  doi:10.1016/j.cpc.2018.10.002, arXiv:1804.01613.
- [KSDT14] V. V. Karasiev, T. Sjostrom, J. Dufty, and S. B. Trickey, "Accurate
  homogeneous electron gas exchange-correlation free energy for local
  spin-density calculations," Phys. Rev. Lett. **112**, 076403 (2014),
  doi:10.1103/PhysRevLett.112.076403.
- [corrKSDT] V. V. Karasiev, J. W. Dufty, and S. B. Trickey, "Nonempirical
  semilocal free-energy density functional for matter under extreme
  conditions," Phys. Rev. Lett. **120**, 076401 (2018),
  doi:10.1103/PhysRevLett.120.076401. The corrected KSDT parameters are in the
  Supplemental Material.
- [KSDT-status] V. V. Karasiev, S. B. Trickey, and J. W. Dufty, "Status of
  free-energy representations for the homogeneous electron gas," Phys. Rev. B
  **99**, 195134 (2019), doi:10.1103/PhysRevB.99.195134. Cross-check for the
  corrected parameters.
- [libxc] Libxc, functional `LDA_XC_CORRKSDT` (ID 318),
  `src/lda_xc_ksdt.c` and `python/lda_exc/lda_xc_ksdt.py`, MPL-2.0,
  https://gitlab.com/libxc/libxc. Open-source implementation of [corrKSDT]; used
  here as the source of the corrected coefficients because the 2018 PRL
  Supplemental Material was not available.
- [TON98] O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura, Phys. Rev. A
  **58**, 2077 (1998), doi:10.1103/PhysRevA.58.2077.
- [BT07] D. A. Batishchev and O. I. Tolstikhin, Phys. Rev. A **75**, 062704
  (2007), doi:10.1103/PhysRevA.75.062704.
- [TON97] O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura, Phys. Rev. Lett.
  **79**, 2026 (1997), doi:10.1103/PhysRevLett.79.2026.

Optional validation references:

- [PZ81] J. P. Perdew and A. Zunger, Phys. Rev. B **23**, 5048 (1981),
  doi:10.1103/PhysRevB.23.5048.
- [PW92] J. P. Perdew and Y. Wang, Phys. Rev. B **45**, 13244 (1992),
  doi:10.1103/PhysRevB.45.13244.
- [Groth17] S. Groth et al., Phys. Rev. Lett. **119**, 135001 (2017),
  doi:10.1103/PhysRevLett.119.135001. Out of scope for implementation; useful
  only for a comparison.

## Units and conventions

- Internal: Hartree and Bohr. Public `average_atom`: `T` in eV, `rho_m` in
  g/cm^3, `A` in u.
- Thermodynamics: internal Hartree; public eV/atom for `F`, `U`; `k_B`/atom for
  `S`; GPa for `P`.
- `k_B = 3.166811563455608e-6` Ha/K; `1 Ha/bohr^3 = 2.9421015696522098e4` GPa.
- `P_e` is the SS23 excess pressure (total minus the ideal-ion `k_BT/V`).
  `P_tot = P_e + k_BT/V`. No ion correlations.
- Negative pressure is returned unchanged, not clamped.

## Interface decisions

### XC model

- Introduce `AbstractXC`.
- Required method: `xc_potential(xc, n, r, T) -> V_xc` (Hartree on the grid).
- Optional method: `xc_free_energy(xc, n, T) -> f_xc` (Hartree per particle on
  the grid), default `nothing`.
- `NoXC` for `xc = nothing`: `V_xc = 0`, `f_xc = 0`.
- `KSDT` implements both with the corrKSDT unpolarized LDA; `xc = :ksdt` is
  sugar for `KSDT()`.
- A raw callable `xc(n, r) -> V_xc` is wrapped as a potential-only XC model; the
  wrapper provides `V_xc` and `f_xc = nothing`.
- `thermodynamics` throws `ArgumentError` when `xc_free_energy === nothing`.
  It must not silently set the XC contributions to zero.
- M3a keeps `xc = nothing` as the default; `:ksdt` is opt-in. Decide whether to
  make `:ksdt` the default after M3c.

### Result and thermodynamics API

- `AverageAtomResult` gains `nuclear_charge::Float64` and `xc::AbstractXC`.
- `thermodynamics(aa::AverageAtomResult)` returns `ThermodynamicProperties`.
- `ThermodynamicProperties` fields:
  - `free_energy`, `internal_energy` (Ha/atom)
  - `entropy` (Ha/K/atom)
  - `electron_pressure`, `ion_ideal_pressure`, `total_pressure` (Ha/bohr^3)
  - `free_energy_electrostatic`, `free_energy_xc`, `internal_energy_kinetic`,
    `entropy_ks`, `entropy_xc`, `pressure_xc`
  - convenience properties: `free_energy_eV`, `internal_energy_eV`,
    `entropy_kB`, `electron_pressure_GPa`, `total_pressure_GPa`
- `F`, `U`, `S` are per atom; `P` is per volume.
- The algebraic identity `F = U - T(S_ks + S_xc)` must hold.

## M3a: corrKSDT XC

Implement the corrKSDT unpolarized LDA:

- `r_s = (3/(4*pi*n))^(1/3)`, `t = T/T_F`,
  `T_F = (3*pi^2*n)^(2/3)/2` in Hartree.
- `f_xc(r_s,t)` is the closed-form KSDT fit: a rational function of
  `sqrt(r_s)` with `t`-dependent coefficients built from rational functions of
  `t^2`, `tanh(1/t)`, `tanh(1/sqrt(t))`, and `exp(-c3/t)`.
- Derivatives (the `t(n)` chain rule is required):
  - `V_xc = f - (r_s/3) f_rs - (2/3) t f_t`
  - `S_xc = -f_t/T_F`
  - `U_xc = f - t f_t`
  - `P_xc = n * (-(r_s/3) f_rs - (2/3) t f_t)`
  - `F_xc = integral n f_xc d^3r`
- Use the corrected parameters from [corrKSDT] (Supplemental Material when
  available), cross-checked against [KSDT-status]. Because the Supplemental
  Material was not available, take the coefficients from [libxc]
  (`LDA_XC_CORRKSDT`). The unpolarized values are:
  ```
  b = (0.342554, 9.141315, 0.448483, 18.553096, 1.05414999729322402165649834296)
  c = (0.875130, -0.256320, 0.953988)
  d = (0.725917, 2.237347, 0.280748, 4.185911, 0.692183)
  e = (0.255415, 0.931933, 0.115398, 17.234117, 0.451437)
  ```
  Keep the source and constants in a comment. If the Supplemental Material
  becomes available, verify the coefficients before changing anything. Do not
  refit.
- No Groth17, no gradient corrections, no spin polarization.
- Stable evaluation at `t -> 0`: `tanh`/`exp` derivatives must not overflow or
  produce `NaN`; use stable forms for `sech^2` and `exp`.
- ForwardDiff is a test-only dependency, used to check the analytic derivatives.

Acceptance criteria:

- [ ] The corrKSDT `f_xc` matches the [libxc] implementation at representative
      `(r_s,t)` points (the 2018 Supplemental Material was unavailable), and the
      analytic high-`T` and `t -> 0` limits hold.
      verify: `just test xc`
- [ ] Analytic `f_rs` and `f_t` match ForwardDiff on
      `r_s in [0.1, 10]`, `t in [1e-6, 10]` (and the `t -> 0` limit), with
      maximum relative error `< 1e-6`.
      verify: `just test xc`
- [ ] `V_xc = f + n df/dn` (including `t(n)`) holds to `1e-6` relative.
      verify: `just test xc`
- [ ] `average_atom(...; xc = :ksdt)` converges at the M2 tolerances on the
      reference case, and `xc = nothing` leaves the M2 result unchanged.
      verify: `just test average_atom xc`

Stop conditions: the [libxc] coefficients cannot be reproduced or differ
materially from [KSDT-status]; the analytic derivatives fail the ForwardDiff
check twice; the XC model needs a new runtime dependency.

## M3b: thermodynamic observables

Implement the Tartarus formulas:

- `F_el = (1/2) integral V_H n d^3r`
- `U_k = integral f(eps,mu) chi(eps) eps d(eps) - integral V_eff n d^3r`
- `S_ks = -integral chi(eps) [f ln f + (1-f) ln(1-f)] d(eps)`
- `F_ks = U_k - T S_ks`
- `F = F_el + F_xc + F_ks`
- `U = F_el + U_xc + U_k`
- `S = S_ks + S_xc`
- `P_e = (2 U_k + F_el)/(3 V) + P_xc`
- `P_ion = k_B T/V`
- `P_tot = P_e + P_ion`
- `P_xc = (1/V) [-F_xc + integral n V_xc d^3r]`

DOS construction:

- Bound states: explicit sum over `eps < E_min` (using `aa.contour.E_min` and a
  `bound_tolerance`, default `1e-8`), with weight `2(2l+1)`.
- Continuum: real-axis integral from `E_min` to a cutoff where
  `f(eps,mu) ~ 1e-10`; panel + 4-point Gauss-Legendre, with extra resolution
  near the Fermi edge.
- `siegert_dos` is effectively per unit energy (the free-electron test compares
  to `V k/pi^2`). Confirm this with the free-electron and neutrality tests and
  fix the docstring if needed.
- Free tail for `l > l_max`: analytic free partial DOS, matching the density
  tail convention (`n_free - sum_{l<=l_max} n_free,l`).
- Start with the Siegert DOS for the continuum. If the `N = 100` criteria fail,
  implement a direct-G DOS from `Im tr G` on the real axis and use it for the
  `:direct_green` backend; `:siegert` keeps the Siegert DOS.

Acceptance criteria:

- [ ] `F = U - T(S_ks + S_xc)` within `1e-8` Ha/atom.
      verify: `just test thermodynamics`
- [ ] DOS neutrality `|integral f chi d(eps) - Z|/Z < 1e-6` at the converged
      `mu` for the reference cases.
      verify: `just test thermodynamics`
- [ ] Free-electron gas limit: `P` within 1%, `S` within `1e-3 k_B/atom`,
      `F,U` within `1e-4` Ha/atom (Sommerfeld included).
      verify: `just test thermodynamics`
- [ ] `N = 50` vs `N = 100`: `P` within 1%, `S` within `1e-3 k_B/atom`,
      `F,U` within `1e-4` Ha/atom. If the Siegert DOS fails, the direct-G DOS
      fallback meets the same criteria.
      verify: `just test thermodynamics`
- [ ] `thermodynamics` with a potential-only XC throws `ArgumentError`.
      verify: `just test thermodynamics`
- [ ] The virial and thermodynamic pressure difference is reported, with no
      hard threshold.
      verify: validation script output.

Stop conditions: the Siegert DOS fails the `N = 100` criteria and the direct-G
DOS is not feasible within this spec; the entropy integral needs a fitted
correction; the identities cannot be met without changing the M2 density.

## M3c: validation

Core unit tests: M3a and M3b criteria.

SS23 physical results (direct-G density + Siegert states; `xc = :ksdt`):

- Fig. 1: Al at `2.7 g/cm^3`, `T = 10/20/30` eV, `l = 0` eigenvalues and DOS.
- Fig. 2: Fe at `8 g/cm^3`, `T = 1` eV, `l = 2` eigenvalues and DOS.
- Fig. 3: Fe effective potential at `8 g/cm^3`, `T = 1` eV, `l = 2`.
- Figs. 4/5: C at `T = 10` eV, selected densities over the reported range.
- Fig. 6: C excess pressure at `T = 10/100/1000` eV, three densities per
  temperature spanning the figure.
- Figs. 7/8: not required. Report the direct-G `N`-independence and the
  `:siegert` `N <= 50` trend as a limitation.

Low-temperature bonus (one case each, proposed conditions):

- Al: `Z = 13`, `A = 27`, `rho_m = 2.7 g/cm^3`, `T = 0.1` eV.
- Na: `Z = 11`, `A = 23`, `rho_m = 0.97 g/cm^3`, `T = 0.1` eV.
- Fe: `Z = 26`, `A = 56`, `rho_m = 7.0 g/cm^3`, `T = 0.2` eV.

These are bonus cases: convergence, positive density, DOS neutrality, and finite
thermodynamic values are required; no literature reference is required.

Runtime budget: the M3c validation suite should finish within about one hour on
the current machine (a faster machine may reduce this). A full paper campaign
over a wide density/temperature grid is a separate run, not part of the test
suite.

Acceptance criteria:

- [ ] SS23 Figs. 1-6 subset agrees qualitatively (DOS peaks, `V_eff` shape,
      eigenvalue tracks) and Fig. 6 quantitatively within a few percent or
      visually within the figure resolution.
      verify: `julia --project test/validation/m3.jl`
- [ ] The low-temperature Al/Na/Fe cases converge and satisfy the DOS
      neutrality/positivity criteria.
      verify: `julia --project test/validation/m3.jl lowt`
- [ ] Formatting and the full test suite pass.
      verify: `treefmt -c && just test`

Stop conditions: the paper figures cannot be digitized and visual comparison is
inconclusive; a validation case needs a condition outside the model
(relativistic, spin, gradients); the runtime budget is exceeded twice.

## Out of scope

- Repairing the Siegert density at `N = 150` (separate M2.6/M4).
- Relativistic effects, spin polarization, gradient corrections (Groth17),
  non-spherical potentials, multi-channel scattering.
- Ion-ion correlations, a full ion equation of state, transport coefficients.
- High-temperature contour entropy optimization (the Siegert DOS is used).
- The Tartarus relativistic lutetium figures.
- A full SS23 figure campaign beyond the M3c subset.

## Files

`src/xc.jl` (new), `src/average_atom.jl`, `src/Siegert.jl`,
`test/test_xc.jl` (new), `test/test_thermodynamics.jl` (new),
`test/test_average_atom.jl`, `test/validation/m3.jl` (new), `README.md`,
`Project.toml` (test-only ForwardDiff), `SPEC.md`, `ISSUE.md`.

## Decisions

- 2026-10-31 -- corrKSDT unpolarized LDA for both `V_xc` and `f_xc`; analytic
  derivatives; ForwardDiff is test-only.
- 2026-10-31 -- The 2018 PRL Supplemental Material was not available. Use the
  corrected corrKSDT coefficients from [libxc] (`LDA_XC_CORRKSDT`), with
  [KSDT-status] as the cross-check. The uncorrected KSDT14 set is not used.
  The low-temperature difference is small (at `r_s=1`, about 0.06% at
  `t=0.01`, 0.8% at `t=0.1`, 1.7% at `t=1`), but corrKSDT is preferred.
- 2026-10-31 -- `AbstractXC` with required `xc_potential` and optional
  `xc_free_energy`; `thermodynamics` errors if `f_xc` is unavailable. A
  potential-only XC model is explicit and cannot silently enter the EOS.
- 2026-10-31 -- `AverageAtomResult` gains `nuclear_charge` and `xc`;
  `thermodynamics(aa)` returns `ThermodynamicProperties`.
- 2026-10-31 -- Units: internal Hartree/Bohr; public eV/atom, `k_B`/atom, GPa.
- 2026-10-31 -- `xc = nothing` remains the default in M3a; `:ksdt` is opt-in.
  Decide on the default after M3c.
- 2026-10-31 -- Thermodynamic DOS: Siegert DOS first, with the direct-G DOS as
  the fallback if the `N = 100` criteria fail.
- 2026-10-31 -- Validation: SS23 Figs. 1-6 subset is the physical target;
  Figs. 7/8 are a documented limitation. Low-temperature Al/Na/Fe are bonus
  cases. The SS23 `N = 150` Siegert-density repair is a separate milestone.
- 2026-10-31 -- `P_e` is the SS23 excess pressure; `P_tot = P_e + k_BT/V`;
  negative pressure is not clamped.

## Confirmed

2026-10-31 -- user: "ありがとう。はい、確定です"
