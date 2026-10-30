# Siegert.jl

Siegert-pseudostate (SPS) tools for finite-range potentials, with a
fixed-potential average-atom (AA) density layer.

The package is non-relativistic and uses Hartree atomic units internally:
\(\varepsilon = k^2/2\), energies and temperatures in Hartree, and radii in
Bohr.

## Implemented physics

### Siegert pseudostates and scattering

- Gauss-Radau-Jacobi DVR basis for angular momentum `l`, following
  Tolstikhin-Ostrovsky-Nakamura (TON98) Appendix C and
  Batishchev-Tolstikhin (BT07) Eqs. (91)-(94).
- Linearized SPS quadratic eigenvalue problem for all `2N + l` states,
  including bound, antibound, resonant, and antiresonant states
  (TON98 Eq. (20)/(C24); BT07 Eq. (31)).
- Raw state type `SPS` and the common `AbstractSPS` interface.
- Scattering matrix from the product formula (BT07 Eq. (85)) and phase
  shifts.

### Normalization and Green's function

- Siegert-state normalization with the surface term
  (SS23 Eq. (6); BT07 Eq. (37)).
- `NormalizedSPS`, which wraps a raw `SPS` and adds normalized coefficients,
  the raw normalization constants, and the DVR basis.
- Partial-wave Green's function
  \(G_l(r,r';k) = \sum_n P_{n,l}(r)P_{n,l}(r')/[k_n(k-k_n)]\)
  (SS23 Eq. (9); BT07 Eq. (73)).
- Density of states from the Green's function
  (SS23 Eqs. (16)-(17)).

### Fixed-potential average-atom density (M1)

- Continuum density from the upper-half-plane contour plus the enclosed
  Matsubara poles (SS23 Eq. (13); S15 Eq. (16) and Appendix B).
- Explicit bound-state sum below `E_min` (default `E_min = 0`).
- Chemical potential from charge neutrality,
  \(Z = 4\pi\int_0^R r^2 n_e(r)\,dr\).
- `fixed_potential_density(V, Z, R, T; N, l_max)`, returning an `AADensity`
  with the density, `μ`, states, common grid, neutrality residual, and contour
  metadata.

M1 does **not** update the effective potential. It is the fixed-potential
input for the self-consistent M2 loop.

### Self-consistent average atom (M2)

- Fixed-point iteration of
  \(V_{\rm eff} = -Z/r + V_H + V_{\rm xc}[n] - V_{\rm xc}[n(R)]\)
  (SS23 Eq. (11)), with the M1 contour + Matsubara density.
- Pluggable exchange-correlation through `xc = nothing` or a callable
  `xc(n, r) -> V_xc[n](r)`; the driver applies the `n(R)` subtraction.
- `average_atom(Z, A, rho_m, T)` accepts `T` in eV, `rho_m` in g/cm^3, and
  `A` in atomic mass units, and returns an `AverageAtomResult` with the
  density, potential, `μ` in eV, neutrality and self-consistency residuals,
  iteration count, convergence flag, states, grid, and contour.
- Optional free-electron subtraction (S15 Eqs. (13)/(14)): each partial wave
  is evaluated as \(n_l - n_l^{\rm free}\) and the analytic uniform free
  density is added back. Validated at fixed potential; the self-consistent
  loop is stable with it only for small `l_max` (see limitations).

## Usage

```julia
using Siegert

# Scattering
V = square_well(-112.5, 1.0)
sps = SPS(19, 0, V, 1.0)
S = scattering_matrix(sps)
delta = 0.5 * atan(imag(S(20.0)), real(S(20.0)))

# Fixed-potential average atom
R = 1.0
T = 0.2                    # Hartree
Veff = r -> r >= R ? 0.0 : -5.0
aa = fixed_potential_density(Veff, 5.0, R, T; N = 30, l_max = 3)

aa.chemical_potential
aa.density
aa.neutrality_residual
aa.contour

# Self-consistent average atom: T in eV, rho_m in g/cm^3, A in u
aa2 = average_atom(5.0, 10.0, 10.0, 0.2; N = 50, l_max = 10)
aa2.converged
aa2.chemical_potential_eV
aa2.self_consistency_residual
```

SPS states are raw by default for scattering-only use. The AA layer wraps them
in `NormalizedSPS`:

```julia
sps = SPS(20, 0, V, R)
nsps = NormalizedSPS(sps)
P = siegert_wavefunction(nsps)
G = siegert_green(nsps, 0.5)
```

## Paper-to-code map

| Reference | Implemented content |
| --- | --- |
| TON98, Phys. Rev. A **58**, 2077 (1998) | DVR basis, SPS algebraic problem, \(l=0\) normalization structure, product formula |
| BT07, Phys. Rev. A **75**, 062704 (2007) | Nonzero-`l` DVR/SPS, normalization Eq. (37), Green's function Eq. (73), product formula Eq. (85) |
| SS23, Phys. Rev. E **107**, 025204 (2023) | Eq. (6) normalization, Eq. (9) Green's function, Eqs. (16)-(17) DOS, Eq. (11) `V_eff`, Eq. (13) contour + Matsubara density |
| S15, High Energy Density Phys. **16**, 18 (2015) | Green's-function AA formulation, Matsubara representation, Eqs. (13)/(14) free-electron subtraction, Appendix B contour construction |

References:

- O. I. Tolstikhin, V. N. Ostrovsky, and H. Nakamura,
  "Siegert pseudostate formulation of scattering theory: One-channel case,"
  *Phys. Rev. A* **58**, 2077-2095 (1998),
  doi:10.1103/PhysRevA.58.2077. Cited as TON98.
- D. A. Batishchev and O. I. Tolstikhin,
  "Siegert pseudostate formulation of scattering theory: Nonzero angular
  momenta in the one-channel case," *Phys. Rev. A* **75**, 062704 (2007),
  doi:10.1103/PhysRevA.75.062704. Cited as BT07.
- C. E. Starrett and N. R. Shaffer, "Average Atom Model with Siegert States,"
  *Phys. Rev. E* **107**, 025204 (2023),
  doi:10.1103/PhysRevE.107.025204, arXiv:2212.06229. Cited as SS23.
- C. E. Starrett, "A Green's function quantum average atom model,"
  *High Energy Density Physics* **16**, 18-22 (2015),
  doi:10.1016/j.hedp.2015.05.001. Cited as S15.

## Current limitations

- The self-consistent loop is practical only for moderate bases (`N` up to
  ~50 with `l_max = 10`, or `N ~ 60` with smaller `l_max`). The M1 contour +
  Matsubara sum uses all `2N + l` Siegert states, and the high-`|k|` states
  have near-null norms whose Float64 cancellation limits the density to
  ~1e-4 relative. Larger `N` can wander or fail to bracket the chemical
  potential.
- The self-consistency threshold is `1e-4` Ha, set by that noise floor rather
  than by a 1e-8 numerical target.
- Free-electron subtraction is validated at fixed potential, but in the
  self-consistent loop the free tail is not accurate enough at `l_max = 10`
  for strongly structured potentials.
- No exchange-correlation model beyond the pluggable interface, and no
  pressure, entropy, free energy, or equation of state.
- Non-relativistic, spherical, single-channel, integer `l`.
- Off-diagonal Green's function \(G_l(r,r';k)\) for \(r\ne r'\) is not exposed.

## Roadmap

### M2: self-consistent average atom

- Implemented: `V_eff` from `-Z/r`, the Hartree potential, and a pluggable
  `V_xc`; fixed-point iteration; `average_atom(Z, A, rho_m, T)` with eV and
  g/cm^3 conversions; optional free-electron subtraction.
- Remaining M1-level work: a better-conditioned density evaluation for larger
  bases (S15-style direct Green function or a Siegert resolvent formulation).

### M3: exchange-correlation, pressure, and equation of state

- Finite-temperature KSDT exchange-correlation.
- Pressure, entropy, free energy, and equation-of-state observables.
- Comparison with the SS23 figures and with Tartarus.

## Development

The project targets Julia 1.13 and provides a Nix dev shell.

```bash
nix develop
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
just test
treefmt -c
```

Tests can be run selectively, for example `just test dvr sps` or
`just test average_atom`.
