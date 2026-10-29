# ISSUE -- Average atom model with Siegert states

## What is missing

`SPS(N, l, V, R)` returns Siegert momenta and DVR coefficients for a
finite-range potential. It does not provide the quantities needed to solve an
average atom (AA) model:

- Siegert states normalized in the sense of Starrett & Shaffer Eq. (6);
- the partial-wave Green's function of Eq. (9);
- the electron density and density of states (DOS) obtained from that Green's
  function (Eqs. (12), (13), (16), (17));
- the chemical potential fixed by charge neutrality;
- the self-consistent effective potential, exchange-correlation, and pressure.

The raw eigenvectors returned by the current linearized eigenproblem are not
normalized for Eq. (6). A direct check on a square well shows that the
orthogonality relation holds with the DVR metric `ρ` plus the surface term, not
with the identity. Without this normalization, the Green's function and every
density built from it carry an unknown state-dependent factor.

## Who feels it, and when

Anyone who wants to use the package for finite-temperature average-atom or
warm-dense-matter electronic structure: bound and continuum states, DOS,
pressure ionization, and equation-of-state work. The current API stops at
scattering observables.

## What was tried

The existing SPS implementation covers non-zero angular momenta, phase shifts,
scattering lengths, and bound-state counts. It was verified against Batishchev
& Tolstikhin. The missing step is not another eigenvalue solver; it is the
normalization and Green's-function layer that turns the SPS spectrum into a
density.

The paper's production method for the continuum density is a contour integral
in the upper half energy plane plus a Matsubara sum. A direct real-axis
continuum integrator would be a different method, would not exercise the
contour code, and is deliberately not used as the reference for the production
path.

## What "fixed" looks like -- in the user's words

"Implement the method of Average Atom Model with Siegert States." (translated)

"Contour + Matsubara should be in the first stage. It is also needed for the
fixed-`V_eff` calculation. I do not want to validate with a different
integration method and introduce an extra source of uncertainty." (translated)

## Open questions

None for the first increment. The scope and decisions are in `SPEC.md`. The
self-consistent AA loop, exchange-correlation, and pressure are later increments
and are out of scope here.
