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

### Fixed-potential average-atom density

- Continuum density from the upper-half-plane contour plus the enclosed
  Matsubara poles (SS23 Eq. (13); S15 Eq. (16) and Appendix B).
- Explicit bound-state sum below `E_min` (default `E_min = 0`).
- Chemical potential from charge neutrality,
  \(Z = 4\pi\int_0^R r^2 n_e(r)\,dr\).
- `fixed_potential_density(V, Z, R, T; N, l_max)`, returning an `AADensity`
  with the density, `μ`, states, common grid, neutrality residual, and contour
  metadata.

The fixed-potential layer does **not** update the effective potential. It is
the input for the self-consistent loop.

### Self-consistent average atom

- Fixed-point iteration of
  \(V_{\rm eff} = -Z/r + V_H + V_{\rm xc}[n] - V_{\rm xc}[n(R)]\)
  (SS23 Eq. (11)).
- The default density backend is the S15-style direct Green function
  (`density_method = :direct_green`): the regular/irregular radial solutions
  are integrated at each contour/Matsubara momentum, so the density does not
  use the near-null high-\(k\) Siegert pseudostates. `density_method =
  :siegert` keeps the contour + Matsubara spectral sum.
- Pluggable exchange-correlation through `xc = nothing` (zero XC), `xc =
  :ksdt` (corrKSDT finite-temperature LDA), an `AbstractXC` model, or a raw
  callable `xc(n, r) -> V_xc[n](r)`; the driver applies the `n(R)`
  subtraction.
- `average_atom(Z, A, rho_m, T)` accepts `T` in eV, `rho_m` in g/cm^3, and
  `A` in atomic mass units, and returns an `AverageAtomResult` with the
  density, potential, `μ` in eV, neutrality and self-consistency residuals,
  iteration count, convergence flag, states, grid, and contour.
- The `:siegert` backend supports free-electron subtraction (S15
  Eqs. (13)/(14)); the direct backend always adds the analytic free tail for
  \(l > l_{\max}\).

### Finite-temperature XC and thermodynamics

- corrKSDT unpolarized finite-temperature LDA for both `V_xc` and the XC free
  energy (Karasiev-Dufty-Trickey, with the corrected coefficients from libxc
  `LDA_XC_CORRKSDT`).
- Analytic `f_rs`, `f_t`, and `V_xc = f + n df/dn`, checked against ForwardDiff
  and libxc reference values.
- `thermodynamics(aa)` returns `ThermodynamicProperties`: free energy,
  internal energy, entropy, and electron/ion/total pressure, decomposed into
  electrostatic, kinetic, and XC terms.
- The electron pressure is the excess pressure (total minus the ideal-ion
  contribution) from the Tartarus continuum/free-electron decomposition:
  bound/core states stay in the density and in `F`/`U`, but they are excluded
  from the pressure virial. The literal all-state virial is a diagnostic
  quantity, not an alternative definition of the excess pressure.
- DOS construction uses bound states below `E_min`, a real-axis continuum
  integral, and the analytic free tail for \(l > l_{\max}\).

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

# Finite-temperature XC and equation of state
aa3 = average_atom(13.0, 27.0, 2.7, 10.0; xc = :ksdt)
th = thermodynamics(aa3)
th.free_energy_eV
th.internal_energy_eV
th.entropy_kB
th.electron_pressure_GPa
th.total_pressure_GPa
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
| S15, High Energy Density Phys. **16**, 18 (2015) | Green's-function AA formulation, Jülich regular/irregular solutions, direct Green function Eqs. (13)/(14), free-electron subtraction, Appendix A matching, Appendix B contour construction |
| KSDT14/corrKSDT, PRL **112**, 076403 (2014); PRL **120**, 076401 (2018) | Finite-temperature corrKSDT unpolarized LDA `f_xc`, `V_xc`, entropy, and pressure |
| Tartarus, Comput. Phys. Commun. **235**, 50 (2019) | Thermodynamic observables, DOS construction, continuum/free-electron pressure decomposition |

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

- The `:siegert` spectral backend is still limited to moderate bases (`N` up
  to ~50 with `l_max = 10`) because the contour + Matsubara sum uses all
  `2N + l` pseudostates and the high-\(|k|\) near-null norms lose Float64
  accuracy. The default `:direct_green` backend avoids this and is stable at
  `N = 100/150`.
- The direct backend's self-consistency threshold is `1e-4` Ha (the physical
  target); the μ solve is a safeguarded secant with `mu_tolerance = 1e-10`
  Ha absolute.
- Low-temperature (`T ~ 0.1 eV`) Al/Na DOS neutrality is limited to about
  `1e-5` relative by near-threshold anti-bound states; the SS23 reference cases
  meet `1e-6`.
- Non-relativistic, spherical, single-channel, integer `l`.
- Off-diagonal Green's function \(G_l(r,r';k)\) for \(r\ne r'\) is not exposed.

## Validation

The SS23 figure comparison and the low-temperature Al/Na/Fe checks live in
`test/validation/ss23.jl`:

```bash
julia --project test/validation/ss23.jl        # SS23 Figs. 1-6 subset
julia --project test/validation/ss23.jl lowt   # low-temperature cases
julia --project test/validation/ss23.jl all    # both
```

The full unit test suite and formatting check are run with `just test` and
`treefmt -c`.

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
