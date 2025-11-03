# Analytic potential functions for electron scattering calculations

"""
    square_well(V₀, a)

Create a finite square well potential function.

Returns a potential V(r) with constant value V₀ inside the well and zero outside:
- V(r) = V₀ for r < a (inside the well)
- V(r) = 0 for r ≥ a (outside the well)

# Arguments
- `V₀`: Potential value inside the well (negative for attractive, atomic units)
- `a`: Well radius (atomic units)

# Returns
A function V(r) that can be evaluated at any radial coordinate r.

# Example
```julia
V = square_well(-112.5, 1.0)  # Attractive well
K̃, ξ, L, H̃, z, Λ, ψ, r = sps_matrices(N, l, a, V)
```
"""
square_well(V₀, a) = r -> r > a ? 0.0 : V₀
