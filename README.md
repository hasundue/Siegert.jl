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
| SS23, Phys. Rev. E **107**, 025204 (2023) | Eq. (6) normalization, Eq. (9) Green's function, Eqs. (16)-(17) DOS, Eq. (13) contour + Matsubara density |
| S15, High Energy Density Phys. **16**, 18 (2015) | Green's-function AA formulation, Matsubara representation, Appendix B contour construction |

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

- Fixed `V_eff` only; no self-consistent AA loop.
- No exchange-correlation model, pressure, entropy, free energy, or equation
  of state.
- No free-electron subtraction / high-`l` acceleration. This currently limits
  the real-axis free-electron DOS test and the accuracy of high angular
  momentum sums.
- Non-relativistic, spherical, single-channel, integer `l`.
- Off-diagonal Green's function \(G_l(r,r';k)\) for \(r\ne r'\) is not exposed.

## Roadmap

### M2: self-consistent average atom

- Build `V_eff` from `-Z/r`, the Hartree potential, and a pluggable `V_xc`.
- Iterate the density, potential, and chemical potential to self-consistency.
- Add free-electron subtraction for the angular-momentum sum.
- Expose `average_atom(Z, A, rho_m, T)`, converting `T` from eV and `rho_m`
  from g/cm^3 at the public API boundary.

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
