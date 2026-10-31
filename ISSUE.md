# ISSUE -- M3: finite-temperature XC and equation-of-state observables

## What's wrong

The M2/M2.5 average-atom (AA) layer computes a self-consistent density and
effective potential, but it stops at the density. It has no finite-temperature
exchange-correlation (XC) model of its own and no thermodynamic observables.
`average_atom` accepts a pluggable `xc` callable that returns only `V_xc`; the
SCF can include an XC potential, but the result carries no `f_xc`, free energy,
entropy, internal energy, or pressure. The M3 roadmap items (finite-temperature
XC, pressure, entropy, free energy, equation of state) are unimplemented.

The gap is not only missing formulas. The SCF and the EOS must use the same XC
free energy; a potential-only XC model is not sufficient for a thermodynamically
consistent pressure and entropy. At low temperature (the liquid-metal regime)
the DOS is sharply structured, and the entropy integral must stay on the real
energy axis or on a contour below the first Matsubara pole. The current code has
no thermodynamic DOS path.

The default density backend is the M2.5 direct-G method, which is stable at
`N = 100/150`. The Siegert spectral density is unreliable above `N ~ 50` because
of near-null normalizations; SS23's `N = 150` Siegert-density convergence figures
(Figs. 7/8) are therefore not reproducible with this package. That is a separate
numerical-method limitation, not something M3 fixes.

## Who feels it, and when

- Anyone who wants an equation of state (pressure, entropy, free energy,
  internal energy) from the AA model.
- Anyone reproducing the physical results of SS23 (DOS, pressure ionization,
  excess pressure).
- The planned low-temperature liquid-metal applications (Al, Na, Fe), where
  finite-temperature XC and a consistent entropy and pressure are required.

## What was tried

- M1: fixed-potential Siegert density with the contour + Matsubara
  representation, DOS, and a neutrality solve. No self-consistency and no XC
  model.
- M2: self-consistent `V_eff`, a pluggable `xc` callable returning only `V_xc`,
  unit conversions, and tests. The XC test verifies only potential-level
  effects; no `f_xc` is available, so no XC free energy, entropy, or pressure
  can be computed.
- M2.5: the S15-style direct-G density backend was adopted as the default and is
  stable at `N = 100/150`; `:siegert` remains for DOS, resonances, and M1 tests.
  The `N = 150` Siegert-density convergence figures of SS23 were diagnosed as
  not reproducible (near-null high-`|k|` states), and the direct-G continuum is
  `N`-independent. This limitation is documented in the M2.5 SPEC.

## What "fixed" looks like -- in the user's words

> "低温（液体金属領域）に適用してみたい"

> "論文の完全な再現というのは言いすぎました" (the physical results are the
> target; the `N`-convergence figures are not)

> "SS23 の N 収束は limitation"

> "Al, Na, Fe 1 ケースずつ欲しい" (low-temperature bonus cases)

The immediate goal is a finite-temperature XC model and a set of thermodynamic
observables that are consistent with the SCF state, reproduce the physical
content of SS23's DOS and pressure results, and are usable at low temperature.
The `N = 150` Siegert-density convergence is a separate milestone.

## Open questions

- Which finite-temperature XC parametrization and which published parameter set
  is authoritative?
- Does the Siegert DOS meet the `N = 100` thermodynamic criteria, or is a
  direct-G DOS fallback required?
- Which SS23 figures are feasible to digitize, and where is only a visual
  comparison possible?
- The exact low-temperature conditions for Al, Na, and Fe.
