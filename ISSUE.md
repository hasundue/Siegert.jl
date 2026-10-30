# ISSUE -- M2.5: stable density evaluation for larger Siegert bases

## What's wrong

M2's self-consistent loop is reliable only for moderate bases: `N = 50` with
`l_max = 10` converges, but `N = 100` and `N = 150` either wander, fail to
bracket the chemical potential, or return an invalid density. The paper results
at `N = 150` cannot be reproduced with the current implementation. The M2 SPEC
records this as a known limitation and stops there.

The cause is **not** a different basis in the papers. SS23 cites TON97 [13],
TON98 [14], and BT07 [15], all of which use the same Jacobi-DVR basis as this
package; SS23's normalization is BT07 Eq. (37), which the package already
implements. The remaining differences from the papers are candidates to test:

- the contour may enclose Siegert poles or approach the top of the
  pseudospectrum (an initial-V diagnostic found one enclosed upper-half pole),
- the high-`|k|` Siegert states have near-null norms and their Float64
  cancellation limits the density (`|norm|` down to 1e-14 at `N = 100`),
- the SPEC's `mixing = 0.4` and `-Z/r + Z/R` initial guess differ from
  Tartarus, which uses `alpha = 0.1` simple mixing (80-100 iterations) and
  Eyert/Anderson acceleration,
- SS23 uses the density-level free-electron acceleration of Tartarus Sec. 3.3
  (arXiv footnote of 2212.06229), which the current `free_electron_subtraction`
  only approximates.

## Who feels it, and when

Anyone who needs a converged average-atom density at a larger basis, or who
tries to reproduce the SS23/Tartarus convergence figures. The current package
silently stops being reliable above `N ~ 50`, while the M2 API accepts `N` as a
free parameter.

## What was tried

- M2 implementation and tests: `just test` passes at `N = 50`, `l_max = 10`.
- `N = 100`, `l_max = 10`, `mixing = 0.1`: bracket failure at
  `mu_guess = 1.97`, `lower = 1.87`.
- `N = 100`, `l_max = 10`, `mixing = 0.2`: no convergence after 250 iterations
  (`residV = 6.3e-3`).
- Near-null norms grow with `N` (73/410 states below 1e-8 at `N = 40`,
  579/1010 at `N = 100`).
- Free-electron subtraction at fixed potential passes the M2 criterion, but the
  self-consistent loop with it converges only for `l_max = 4, 6`, not 10.
- The V=0 free Siegert sum is unusable as a reference (1e11 scale); the free
  reference must be analytic (Bessel/Hankel), as in S15 Eq. (14).

## What "fixed" looks like -- in the user's words

> "まず 1. で一区切りして、4. -> 2. -> 3. の順で行う M2.5 を起票"

> "M1.5 には消極的です" (the basis is not the difference)

> "微妙な数値で通らない場合は無理せず中断して相談してください"

## Open questions

- Does the contour enclose Siegert poles for a converged M2 potential, and does
  an enclosed pole explain the `N = 100` garbage?
- Is the dominant error contour-related, near-null cancellation, or mixing?
- Does a faithful Tartarus/SS23 density-level free correction stabilize
  `N = 100/150`?
- If needed, is `alpha = 0.1` or Eyert/Anderson acceptable as a SPEC change?
