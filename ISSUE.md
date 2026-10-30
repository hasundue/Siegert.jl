# ISSUE -- Self-consistent average atom from Siegert states

## What's wrong

M1 can compute the electron density and chemical potential for a fixed
effective potential `V_eff(r)`, but it cannot update `V_eff` from the density.
It is therefore not yet an average-atom model: the Hartree potential and any
exchange-correlation potential are absent, the density is not self-consistent,
and there is no public API that accepts a material density and temperature.

Users must currently iterate outside the package, construct `V_eff` by hand,
and cannot obtain an equilibrium average-atom density or chemical potential.

## Who feels it, and when

Anyone who wants to use the package for finite-temperature average-atom or
warm-dense-matter electronic structure. M1 is only an input layer. Without a
self-consistent loop, the package cannot answer the basic average-atom question:
what density and chemical potential correspond to a given `Z`, mass density,
and temperature?

## What was tried

M1 implemented the Siegert normalization, Green's function, DOS, contour +
Matsubara density, and the fixed-potential driver. Those tests pass. The M1
driver deliberately does not update `V_eff`, and the M1 SPEC left self-
consistency, exchange-correlation, and pressure to later increments.

The roadmap in the archived M1 SPEC lists M2 as a self-consistent AA loop with
a pluggable `V_xc`, free-electron subtraction, and a public
`average_atom(Z, A, rho_m, T)`. That roadmap has no acceptance criteria,
unit decisions, convergence criteria, or stop conditions yet.

## What "fixed" looks like -- in the user's words

> "M2 の keel を"

> "XC を pluggable にするのはとても良いと思う"

> "返り値は eV にしておきましょうか。必要ならあとで単位周りのインタフェースを整備する感じで。"

> "微妙な数値で通らない場合は無理せず中断して相談してください。"

## Open questions

None for the spec. Decisions and acceptance criteria are in `SPEC.md`.
Free-electron subtraction is a separate M2 sub-increment; if its free reference
cannot be defined consistently, that is a blocker rather than a reason to fit a
factor.
