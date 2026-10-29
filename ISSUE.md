# ISSUE -- Non-zero angular momentum support

## What's wrong

`SPS(N, l, ...)` is only correct for `l = 0`. Angular momentum already appears in
the DVR basis and the centrifugal potential, but the parts that make the method
work are still s-wave specific:

- `src/sps.jl` fixes the boundary parameter to `b = -1`.
- `src/sps.jl` linearizes the quadratic eigenproblem to dimension `2N`.
- `src/spectra.jl` applies the product formula to the `2N` momenta.

The s-wave boundary condition is not the outgoing-wave boundary condition for
`l > 0`, and the non-zero-`l` problem has `2N + l` solutions. Using the current
code with `l = 1, 2, ...` therefore gives values that should not be trusted.

## Who feels it, and when

Anyone using the package for `p`-, `d`- or higher partial waves: phase shifts,
the scattering matrix, and Siegert poles. The limitation blocks the one-channel
radial part of any three-dimensional scattering calculation.

## What was tried

The current code has partial `l` plumbing from the original Tolstikhin,
Ostrovsky, and Nakamura (TON98) formulation, but TON98 explicitly leaves
non-zero angular momenta as an open problem. Batishchev and Tolstikhin (2007)
give the generalization for non-zero angular momenta in the one-channel case.
No implementation of that generalization exists in this repository yet.

## What "fixed" looks like -- in the user's words

"SPS with non-zero angular momenta should produce the correct Siegert states and
phase shifts; the existing s-wave results must continue to pass."

## Open questions

None. The interview decisions are recorded in `SPEC.md`.
