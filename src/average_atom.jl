# Average-atom density from Siegert states.
#
# References
# - C. E. Starrett and N. R. Shaffer, "Average Atom Model with Siegert States,"
#   Phys. Rev. E 107, 025204 (2023). doi:10.1103/PhysRevE.107.025204 [SS23]
# - C. E. Starrett, "A Green's function quantum average atom model,"
#   High Energy Density Physics 16, 18-22 (2015).
#   doi:10.1016/j.hedp.2015.05.001 [S15]
# - D. A. Batishchev and O. I. Tolstikhin, "Siegert pseudostate formulation of
#   scattering theory: Nonzero angular momenta in the one-channel case,"
#   Phys. Rev. A 75, 062704 (2007). doi:10.1103/PhysRevA.75.062704 [BT07]

# Physical constants and unit conversions for the public API. Internal
# calculations use Hartree and Bohr; the public `average_atom` accepts eV,
# g/cm^3, and atomic mass units.
const HARTREE_TO_EV = 27.211386245988
const ATOMIC_MASS_UNIT_KG = 1.66053906660e-27
const BOHR_RADIUS_M = 5.29177210903e-11
const GRAM_PER_CM3_TO_KG_PER_M3 = 1000.0


"""
Parameters of the upper-half-plane energy contour used by SS23 Eq. (13).

Fields:
- `E_min`: Lower end of the contour on the real energy axis
- `E_max`: Upper end of the contour on the real energy axis
- `E_I`: Imaginary part of the horizontal section
- `N_mat`: Number of Matsubara poles enclosed by the contour
- `energies`: Enclosed Matsubara energies `μ + iπ(2j-1)T`
"""
struct Contour
    E_min::Float64
    E_max::Float64
    E_I::Float64
    N_mat::Int
    energies::Vector{ComplexF64}
end

"""
Normalized Siegert state.

Wraps a raw `SPS` and adds the normalized coefficients, raw normalization
constants, and the DVR basis. The average-atom layer uses this type.

Fields:
- `sps`: Raw eigenproblem result
- `coefficients`: Normalized DVR coefficients
- `normalizations`: Raw normalization constants `N_n`
- `dvr`: DVR basis used to build and evaluate the state
"""
struct NormalizedSPS <: AbstractSPS
    sps::SPS
    coefficients::Matrix{ComplexF64}
    normalizations::Vector{ComplexF64}
    dvr::DVR
end

# Accessors shared with the raw SPS interface.
momenta(nsps::NormalizedSPS) = momenta(nsps.sps)
coefficients(nsps::NormalizedSPS) = nsps.coefficients
grid(nsps::NormalizedSPS) = grid(nsps.sps)
box_size(nsps::NormalizedSPS) = box_size(nsps.sps)
angular_momentum(nsps::NormalizedSPS) = angular_momentum(nsps.sps)
hamiltonian(nsps::NormalizedSPS) = hamiltonian(nsps.sps)
metric(nsps::NormalizedSPS) = metric(nsps.sps)
boundary(nsps::NormalizedSPS) = boundary(nsps.sps)
dvr(nsps::NormalizedSPS) = nsps.dvr

"""
Compute the SPS normalization constants and rescale the DVR coefficients.

The normalization follows SS23 Eq. (6) and BT07 Eq. (37). The radial
wavefunction is `P_n(r) = r * Σ_i c_{i,n} π_i(x)`, so the volume term is
`(a/2) c_n^T ρ c_n` and the surface value is `P_n(a) = a π(a)^T c_n`.

Returns the raw normalization constants `N_n`.
"""
function _normalize_coefficients!(
    coefficients::Matrix{ComplexF64},
    ks::Vector{ComplexF64},
    rho::Matrix{Float64},
    basis::DVR,
    a::Real,
    l::Integer,
)
    z = reverse_bessel_zeros(l)
    π_at_boundary = basis.basis[:, end]
    norm = Vector{ComplexF64}(undef, length(ks))
    factors = Vector{ComplexF64}(undef, length(ks))

    for i in eachindex(ks)
        c = view(coefficients, :, i)
        k = ks[i]
        P = a * sum(π_at_boundary .* c)
        surface = one(ComplexF64)
        for zp in z
            surface += zp / (im * k * a + zp)^2
        end
        interior = (a / 2) * sum(c .* (rho * c))
        norm[i] = interior + im * P^2 / (2 * k) * surface
        factors[i] = sqrt(norm[i])
    end

    coefficients ./= reshape(factors, 1, :)
    return norm
end

"""
    NormalizedSPS(sps::SPS) -> NormalizedSPS
    NormalizedSPS(sps::SPS, dvr::DVR) -> NormalizedSPS

Wrap a raw `SPS` and apply SS23 Eq. (6) / BT07 Eq. (37) normalization. The
optional `dvr` avoids rebuilding the basis when it is already available.
"""
function NormalizedSPS(sps::SPS, basis::DVR)
    N = size(sps.coefficients, 1)
    size(basis.basis, 1) == N ||
        throw(DimensionMismatch("DVR basis does not match the SPS dimension"))
    coefficients = copy(sps.coefficients)
    norm = _normalize_coefficients!(
        coefficients,
        sps.momenta,
        sps.metric,
        basis,
        sps.box_size,
        sps.angular_momentum,
    )
    return NormalizedSPS(sps, coefficients, norm, basis)
end

function NormalizedSPS(sps::SPS)
    basis = DVR(size(sps.coefficients, 1), sps.angular_momentum)
    return NormalizedSPS(sps, basis)
end

"""
Result of a fixed-potential average-atom density calculation.

Fields:
- `density`: Electron density `n_e(r)` on `grid`
- `chemical_potential`: Chemical potential `μ` in Hartree
- `states`: Normalized Siegert states for each angular momentum
- `grid`: Common radial grid
- `radius`: Ion-sphere radius `R`
- `temperature`: Temperature in Hartree
- `neutrality_residual`: `4π ∫ r² n_e dr - Z` at the returned `μ`
- `contour`: Contour parameters used for the final density
"""
struct AADensity
    density::Vector{Float64}
    chemical_potential::Float64
    states::Vector{NormalizedSPS}
    grid::Vector{Float64}
    radius::Float64
    temperature::Float64
    neutrality_residual::Float64
    contour::Contour
end

"""
Return the `l = 0` Gauss-Radau grid and the weights for `∫_0^R f(r) dr`.

The same grid is used for every angular momentum because the DVR grids depend
on `l`. For a `NormalizedSPS` with `l = 0` its own grid is the common grid.
"""
function common_radial_grid(nsps::NormalizedSPS)
    N = size(coefficients(nsps), 1)
    basis = angular_momentum(nsps) == 0 ? dvr(nsps) : DVR(N, 0)
    r = x_to_r.(basis.grid, box_size(nsps))
    weights = (box_size(nsps) / 2) .* basis.weights
    return r, weights
end

function common_radial_grid(states::AbstractVector{<:NormalizedSPS})
    isempty(states) && throw(ArgumentError("states must not be empty"))
    return common_radial_grid(first(states))
end

"""
Return the weights for the volume integral `∫_0^R r² f(r) dr`.
"""
volume_weights(r::AbstractVector{<:Real}, weights::AbstractVector{<:Real}) =
    (r .^ 2) .* weights

"""
    siegert_wavefunction(nsps::NormalizedSPS) -> Matrix{ComplexF64}

Evaluate all normalized radial Siegert states of `nsps` on `grid(nsps)`.

The returned matrix has one row per grid point and one column per state, with
`P_{n,l}(r) = r Σ_i c_{i,n} π_i(x(r))`.
"""
function siegert_wavefunction(nsps::NormalizedSPS)
    return siegert_wavefunction(nsps, grid(nsps))
end

siegert_wavefunction(sps::SPS) =
    throw(ArgumentError("wrap the raw SPS in NormalizedSPS first"))
siegert_wavefunction(sps::SPS, r::AbstractVector) =
    throw(ArgumentError("wrap the raw SPS in NormalizedSPS first"))

"""
    siegert_wavefunction(nsps::NormalizedSPS, r::AbstractVector) -> Matrix{ComplexF64}

Evaluate all normalized radial Siegert states of `nsps` on the radii `r`.

`r` must lie in `[0, R]`. The returned matrix has one row per radius and one
column per state.
"""
function siegert_wavefunction(nsps::NormalizedSPS, r::AbstractVector)
    all(0 .<= r .<= box_size(nsps)) || throw(DomainError(r, "radii must lie in [0, R]"))
    l = angular_momentum(nsps)
    C = coefficients(nsps)
    N = size(C, 1)
    x = 2 .* r ./ box_size(nsps) .- 1

    w = jacobi_weight.(x, 0, 2l)
    h = jacobi_normsq.(0:(N-1), 0, 2l)
    P = jacobi_polynomials(x, 0, 2l, N)
    φ = sqrt.(reshape(w, 1, :) ./ reshape(h, :, 1)) .* P
    π = transpose(dvr(nsps).transformation) * φ
    values = transpose(C) * π
    return ComplexF64.(transpose(reshape(r, 1, :) .* values))
end

"""
    siegert_green(nsps::NormalizedSPS, k::Number; r=grid(nsps)) -> Vector{ComplexF64}

Diagonal partial-wave Green's function
`G_l(r,r;k) = Σ_n P_{n,l}(r)^2 / [k_n (k - k_n)]` on the radii `r`.
"""
function siegert_green(nsps::NormalizedSPS, k::Number; r = grid(nsps))
    P = siegert_wavefunction(nsps, r)
    return _green_diagonal(P, momenta(nsps), k)
end

function _green_diagonal(P::AbstractMatrix, ks::AbstractVector, k::Number)
    state_weights = @. 1 / (ks * (k - ks))
    return (P .^ 2) * state_weights
end

"""
    siegert_dos(nsps::NormalizedSPS, k::Number; r, weights) -> Float64

Partial density of states `χ_l(k) = -(1/π) Im ∫_0^R dr 2(2l+1) G_l(r,r;k)`.

When `r` and `weights` are omitted, the common `l = 0` grid is used. The
weights must integrate over `dr`, not over `r² dr`.
"""
function siegert_dos(
    nsps::NormalizedSPS,
    k::Number;
    r::Union{Nothing,AbstractVector} = nothing,
    weights::Union{Nothing,AbstractVector} = nothing,
)
    if r === nothing
        r, default_weights = common_radial_grid(nsps)
        weights === nothing && (weights = default_weights)
    end
    weights === nothing && throw(ArgumentError("weights must be given together with r"))
    length(r) == length(weights) ||
        throw(DimensionMismatch("r and weights must have the same length"))
    G = siegert_green(nsps, k; r = r)
    return real(-(1 / π) * imag(2 * (2angular_momentum(nsps) + 1) * sum(weights .* G)))
end

"""
    siegert_dos(states::AbstractVector{<:NormalizedSPS}, k::Number; r, weights)
    -> Float64

Total density of states, summed over the given partial waves.
"""
function siegert_dos(
    states::AbstractVector{<:NormalizedSPS},
    k::Number;
    r::Union{Nothing,AbstractVector} = nothing,
    weights::Union{Nothing,AbstractVector} = nothing,
)
    if r === nothing
        r, default_weights = common_radial_grid(states)
        weights === nothing && (weights = default_weights)
    end
    weights === nothing && throw(ArgumentError("weights must be given together with r"))
    return sum(siegert_dos(nsps, k; r = r, weights = weights) for nsps in states)
end

"""
    fermi_dirac(ε, μ, T)

Fermi-Dirac occupation factor `1 / (1 + exp((ε - μ)/T))`, evaluated in a form
that avoids overflow. Complex `ε` is supported for contour integration.
"""
function fermi_dirac(ε::Real, μ::Real, T::Real)
    T > 0 || throw(ArgumentError("T must be positive"))
    x = (ε - μ) / T
    if x > 0
        e = exp(-x)
        return e / (1 + e)
    else
        return 1 / (1 + exp(x))
    end
end

function fermi_dirac(ε::Complex, μ::Real, T::Real)
    T > 0 || throw(ArgumentError("T must be positive"))
    x = (ε - μ) / T
    if real(x) > 0
        e = exp(-x)
        return e / (1 + e)
    else
        return 1 / (1 + exp(x))
    end
end

"""
    matsubara_energies(μ, T; E_min=-Inf, E_max=Inf, E_I=1.0) -> Vector{ComplexF64}

Matsubara energies `z_j = μ + iπ(2j-1)T` that lie inside the contour limits.
"""
function matsubara_energies(
    μ::Real,
    T::Real;
    E_min::Real = -Inf,
    E_max::Real = Inf,
    E_I::Real = 1.0,
)
    T > 0 || throw(ArgumentError("T must be positive"))
    (E_min < μ < E_max) || return ComplexF64[]
    E_I > 0 || throw(ArgumentError("E_I must be positive"))

    energies = ComplexF64[]
    j = 1
    while π * (2j - 1) * T < E_I
        push!(energies, μ + im * π * (2j - 1) * T)
        j += 1
    end
    return energies
end

"""
    contour_parameters(μ, T; E_min=0.0, E_max=nothing, E_I=nothing) -> Contour

Construct the S15 Appendix B contour parameters for the given `μ` and `T`.
"""
function contour_parameters(
    μ::Real,
    T::Real;
    E_min::Real = 0.0,
    E_max::Union{Nothing,Real} = nothing,
    E_I::Union{Nothing,Real} = nothing,
)
    T > 0 || throw(ArgumentError("T must be positive"))
    E_min = Float64(E_min)

    if E_max === nothing
        E_max = max(0.0, μ + T * log(1e10 - 1))
    end
    E_max = Float64(E_max)
    E_max > E_min || return Contour(E_min, E_max, 1.0, 0, ComplexF64[])

    if E_I === nothing
        target = 1.0
        if !(E_min < μ < E_max) || π * T >= target
            E_I = target
        else
            n = max(1, round(Int, target / (2π * T)))
            E_I = 2n * π * T
        end
    else
        E_I = Float64(E_I)
        E_I > 0 || throw(ArgumentError("E_I must be positive"))
    end

    energies = matsubara_energies(μ, T; E_min = E_min, E_max = E_max, E_I = E_I)
    return Contour(E_min, E_max, Float64(E_I), length(energies), energies)
end

function gauss_legendre4()
    nodes =
        (-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
    weights =
        (0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538)
    return nodes, weights
end

function _append_legendre_segment!(points, measures, z0, z1, panels)
    nodes, weights = gauss_legendre4()
    for panel = 1:panels
        center = -1 + (2panel - 1) / panels
        half = 1 / panels
        for q = 1:4
            t = center + half * nodes[q]
            push!(points, z0 + (z1 - z0) * (t + 1) / 2)
            push!(measures, (z1 - z0) / 2 * half * weights[q])
        end
    end
    return points, measures
end

"""
    contour_points(contour; n_vert=1, n_horiz=20, n_return=20)

Return the three-sided S15 contour points and the complex energy measures `dε`.

The real-axis segment is not part of the contour; it is the original
integration path that the contour replaces.
"""
function contour_points(
    contour::Contour;
    n_vert::Integer = 1,
    n_horiz::Integer = 20,
    n_return::Integer = 20,
)
    n_vert > 0 || throw(ArgumentError("n_vert must be positive"))
    n_horiz > 0 || throw(ArgumentError("n_horiz must be positive"))
    n_return > 0 || throw(ArgumentError("n_return must be positive"))
    points = ComplexF64[]
    measures = ComplexF64[]
    contour.E_max > contour.E_min || return points, measures

    points, measures = _append_legendre_segment!(
        points,
        measures,
        contour.E_max,
        contour.E_max + im * contour.E_I,
        n_vert,
    )
    points, measures = _append_legendre_segment!(
        points,
        measures,
        contour.E_max + im * contour.E_I,
        contour.E_min + im * contour.E_I,
        n_horiz,
    )
    points, measures = _append_legendre_segment!(
        points,
        measures,
        contour.E_min + im * contour.E_I,
        contour.E_min,
        n_return,
    )
    return points, measures
end

function _state_data(states::AbstractVector{<:NormalizedSPS}, r::AbstractVector)
    return [
        (momenta(nsps), siegert_wavefunction(nsps, r) .^ 2, angular_momentum(nsps)) for
        nsps in states
    ]
end

function _is_bound_state(k::ComplexF64, tolerance::Real)
    return abs(real(k)) <= tolerance && imag(k) > 0
end

function _add_bound_density!(density, data, r2, μ, T, E_min, bound_tolerance)
    factor = 1 / (4π)
    for (ks, P2, l) in data
        for i in eachindex(ks)
            _is_bound_state(ks[i], bound_tolerance) || continue
            ε = real(ks[i]^2) / 2
            ε < E_min || continue
            occupation = fermi_dirac(ε, μ, T)
            density .+= real.(occupation .* (2 * (2l + 1) * factor) .* (P2[:, i] ./ r2))
        end
    end
    return density
end

function _add_contour_density!(density, data, r2, contour, μ, T; n_vert, n_horiz, n_return)
    points, measures =
        contour_points(contour; n_vert = n_vert, n_horiz = n_horiz, n_return = n_return)
    isempty(points) && return density

    factor = 1 / (4π)
    accumulator = zeros(ComplexF64, length(density))
    for (ε, dε) in zip(points, measures)
        k = sqrt(2 * ε)
        occupation = fermi_dirac(ε, μ, T)
        for (ks, P2, l) in data
            G = P2 * (1 ./ (ks .* (k .- ks)))
            accumulator .+= occupation * dε * (2 * (2l + 1) * factor) .* (G ./ r2)
        end
    end
    density .+= (1 / π) .* imag.(accumulator)
    return density
end

function _add_matsubara_density!(density, data, r2, contour, T)
    factor = 1 / (4π)
    for z in contour.energies
        k = sqrt(2 * z)
        for (ks, P2, l) in data
            G = P2 * (1 ./ (ks .* (k .- ks)))
            density .+= 2 * T * real.((2 * (2l + 1) * factor) .* (G ./ r2))
        end
    end
    return density
end

# Spherical Bessel and outgoing spherical Hankel functions, used by the
# analytic free-electron reference of S15 Eq. (14).
_spherical_j(l::Integer, z) = sqrt.(π ./ (2 .* z)) .* besselj.(l + 1 / 2, z)
_spherical_h1(l::Integer, z) = sqrt.(π ./ (2 .* z)) .* besselh.(l + 1 / 2, 1, z)

"""
    free_electron_density(μ, T) -> Float64

Uniform free-electron density `(1/π²) ∫_0^∞ k² f(k²/2) dk` in bohr^-3.
"""
function free_electron_density(μ::Real, T::Real)
    T > 0 || throw(ArgumentError("T must be positive"))
    k_max = sqrt(2 * (max(0.0, Float64(μ)) + 60T))
    nodes, weights = gauss_legendre4()
    panels = 200
    integral = 0.0
    for panel = 1:panels
        a = (panel - 1) / panels * k_max
        b = panel / panels * k_max
        half = (b - a) / 2
        center = (a + b) / 2
        for q = 1:4
            k = center + half * nodes[q]
            integral += half * weights[q] * k^2 * fermi_dirac(k^2 / 2, μ, T)
        end
    end
    return integral / π^2
end

# Partial free-electron density from the analytic free Green function
# `G_l^0 = -2 i k r² j_l(kr) h^+_l(kr)` (S15 Eq. 14). The factor 2 comes from
# the spin degeneracy already present in `2(2l + 1)`.
function _free_partial_density(r, μ, T, l, contour, n_vert, n_horiz, n_return)
    points, measures = contour_points(contour; n_vert, n_horiz, n_return)
    factor = 2 * (2l + 1) / (4π)
    density = zeros(Float64, length(r))
    accumulator = zeros(ComplexF64, length(r))
    for (ε, dε) in zip(points, measures)
        k = sqrt(2ε)
        G = -2im * k .* _spherical_j(l, k .* r) .* _spherical_h1(l, k .* r)
        accumulator .+= fermi_dirac(ε, μ, T) * dε * factor .* G
    end
    density .+= (1 / π) .* imag.(accumulator)
    for z in contour.energies
        k = sqrt(2z)
        G = -2im * k .* _spherical_j(l, k .* r) .* _spherical_h1(l, k .* r)
        density .+= 2T * real.(factor .* G)
    end
    return density
end

function _density_and_contour(
    data,
    r::AbstractVector,
    μ::Real,
    T::Real;
    E_min::Real = 0.0,
    E_max::Union{Nothing,Real} = nothing,
    E_I::Union{Nothing,Real} = nothing,
    n_vert::Integer = 1,
    n_horiz::Integer = 20,
    n_return::Integer = 20,
    bound_tolerance::Real = 1e-8,
    free_electron_subtraction::Bool = false,
)
    all(>(0), r) || throw(DomainError(r, "radii must be positive"))
    contour = contour_parameters(μ, T; E_min = E_min, E_max = E_max, E_I = E_I)
    r2 = r .^ 2
    density = zeros(Float64, length(r))

    _add_bound_density!(density, data, r2, μ, T, E_min, bound_tolerance)
    _add_contour_density!(
        density,
        data,
        r2,
        contour,
        μ,
        T;
        n_vert = n_vert,
        n_horiz = n_horiz,
        n_return = n_return,
    )
    _add_matsubara_density!(density, data, r2, contour, T)
    if free_electron_subtraction
        l_values = sort(unique(l for (_, _, l) in data))
        density .+= free_electron_density(μ, T)
        for l in l_values
            density .-=
                _free_partial_density(r, μ, T, l, contour, n_vert, n_horiz, n_return)
        end
    end
    return density, contour
end

"""
    siegert_density(states, μ, T; E_min=0, r, weights, ...) -> Vector{Float64}

Electron density `n_e(r)` from SS23 Eq. (13) on the common radial grid.

Bound states with energy below `E_min` are added explicitly; the continuum is
evaluated with the three-sided contour plus the enclosed Matsubara poles.
"""
function siegert_density(
    states::AbstractVector{<:NormalizedSPS},
    μ::Real,
    T::Real;
    E_min::Real = 0.0,
    E_max::Union{Nothing,Real} = nothing,
    E_I::Union{Nothing,Real} = nothing,
    r::Union{Nothing,AbstractVector} = nothing,
    weights::Union{Nothing,AbstractVector} = nothing,
    wavefunctions::Union{Nothing,AbstractVector} = nothing,
    n_vert::Integer = 1,
    n_horiz::Integer = 20,
    n_return::Integer = 20,
    bound_tolerance::Real = 1e-8,
    free_electron_subtraction::Bool = false,
)
    if r === nothing
        r, default_weights = common_radial_grid(states)
        weights === nothing && (weights = default_weights)
    end
    weights === nothing && throw(ArgumentError("weights must be given together with r"))
    length(r) == length(weights) ||
        throw(DimensionMismatch("r and weights must have the same length"))

    if wavefunctions === nothing
        data = _state_data(states, r)
    else
        length(wavefunctions) == length(states) ||
            throw(DimensionMismatch("wavefunctions must match states"))
        data = [
            (momenta(states[i]), wavefunctions[i] .^ 2, angular_momentum(states[i])) for
            i in eachindex(states)
        ]
    end

    density, _ = _density_and_contour(
        data,
        r,
        μ,
        T;
        E_min = E_min,
        E_max = E_max,
        E_I = E_I,
        n_vert = n_vert,
        n_horiz = n_horiz,
        n_return = n_return,
        bound_tolerance = bound_tolerance,
        free_electron_subtraction = free_electron_subtraction,
    )
    return density
end

"""
    chemical_potential(states, Z, R, T; E_min=0, ...) -> (μ, density, residual)

Solve `Z = 4π ∫_0^R r² n_e(r) dr` for the chemical potential by bisection with
an expanding bracket. The residual is `4π ∫ r² n_e dr - Z`.
"""
function chemical_potential(
    states::AbstractVector{<:NormalizedSPS},
    Z::Real,
    R::Real,
    T::Real;
    E_min::Real = 0.0,
    E_max::Union{Nothing,Real} = nothing,
    E_I::Union{Nothing,Real} = nothing,
    n_vert::Integer = 1,
    n_horiz::Integer = 20,
    n_return::Integer = 20,
    bound_tolerance::Real = 1e-8,
    mu_tolerance::Real = 1e-10,
    max_iterations::Integer = 100,
    mu_guess::Union{Nothing,Real} = nothing,
    free_electron_subtraction::Bool = false,
)
    isempty(states) && throw(ArgumentError("states must not be empty"))
    Z > 0 || throw(ArgumentError("Z must be positive"))
    R > 0 || throw(ArgumentError("R must be positive"))
    T > 0 || throw(ArgumentError("T must be positive"))

    r, quadrature_weights = common_radial_grid(states)
    volume = volume_weights(r, quadrature_weights)
    data = _state_data(states, r)

    function evaluate(mu)
        density, contour = _density_and_contour(
            data,
            r,
            mu,
            T;
            E_min = E_min,
            E_max = E_max,
            E_I = E_I,
            n_vert = n_vert,
            n_horiz = n_horiz,
            n_return = n_return,
            bound_tolerance = bound_tolerance,
            free_electron_subtraction = free_electron_subtraction,
        )
        electron_number = 4π * sum(volume .* density)
        return electron_number - Z, density, contour
    end

    energies = reduce(vcat, (real.(momenta(nsps) .^ 2 ./ 2) for nsps in states))
    if mu_guess === nothing
        lower = minimum(energies) - 10T
        upper = maximum(energies) + 10T
    else
        # Warm start: keep the bracket near the previous chemical potential.
        # The full spectrum's maximum energy can be far above the physical μ;
        # evaluating the contour there may reach the top of the pseudospectrum
        # and return garbage.
        width = max(0.1, 5T, 0.01 * abs(mu_guess))
        lower = Float64(mu_guess) - width
        upper = Float64(mu_guess) + width
    end
    f_lower, _, _ = evaluate(lower)
    f_upper, _, _ = evaluate(upper)

    scale = max(1.0, upper - lower)
    iterations = 0
    while (!isfinite(f_lower) || f_lower > 0) && iterations < max_iterations
        lower -= scale
        scale *= 2
        f_lower, _, _ = evaluate(lower)
        iterations += 1
    end
    scale = max(1.0, upper - lower)
    while (!isfinite(f_upper) || f_upper < 0) && iterations < max_iterations
        upper += scale
        scale *= 2
        f_upper, _, _ = evaluate(upper)
        iterations += 1
    end
    (isfinite(f_lower) && isfinite(f_upper) && f_lower <= 0 && f_upper >= 0) || throw(
        ArgumentError(
            "could not bracket the chemical potential (mu_guess = $mu_guess, " *
            "lower = $lower, f_lower = $f_lower, upper = $upper, f_upper = $f_upper)",
        ),
    )

    for _ = 1:max_iterations
        midpoint = (lower + upper) / 2
        f_mid, _, _ = evaluate(midpoint)
        if f_mid > 0
            upper = midpoint
        else
            lower = midpoint
        end
        (upper - lower) < mu_tolerance && break
    end

    mu = (lower + upper) / 2
    residual, density, _ = evaluate(mu)
    return mu, density, residual
end

"""
    fixed_potential_density(V, Z, R, T; N=150, l_max=10, ...) -> AADensity

Compute the fixed-potential M1 average-atom density for `V`, nuclear charge
`Z`, ion-sphere radius `R`, and temperature `T` in Hartree.

With `free_electron_subtraction = true` each partial wave is evaluated as
`n_l - n_l^free` and the analytic free density is added back (S15
Eqs. (13)/(14)); the total density is unchanged but the `l` sum converges
faster.
"""
function fixed_potential_density(
    V::Function,
    Z::Real,
    R::Real,
    T::Real;
    N::Integer = 150,
    l_max::Integer = 10,
    E_min::Real = 0.0,
    E_max::Union{Nothing,Real} = nothing,
    E_I::Union{Nothing,Real} = nothing,
    n_vert::Integer = 1,
    n_horiz::Integer = 20,
    n_return::Integer = 20,
    bound_tolerance::Real = 1e-8,
    mu_tolerance::Real = 1e-10,
    max_iterations::Integer = 100,
    mu_guess::Union{Nothing,Real} = nothing,
    free_electron_subtraction::Bool = false,
)
    N > 0 || throw(ArgumentError("N must be positive"))
    l_max >= 0 || throw(ArgumentError("l_max must be non-negative"))
    states = map(0:l_max) do l
        sps, basis = _solve_sps(N, l, V, R)
        NormalizedSPS(sps, basis)
    end
    mu, density, residual = chemical_potential(
        states,
        Z,
        R,
        T;
        E_min = E_min,
        E_max = E_max,
        E_I = E_I,
        n_vert = n_vert,
        n_horiz = n_horiz,
        n_return = n_return,
        bound_tolerance = bound_tolerance,
        mu_tolerance = mu_tolerance,
        max_iterations = max_iterations,
        mu_guess = mu_guess,
        free_electron_subtraction = free_electron_subtraction,
    )
    grid, _ = common_radial_grid(states)
    contour = contour_parameters(mu, T; E_min = E_min, E_max = E_max, E_I = E_I)
    return AADensity(density, mu, states, grid, Float64(R), Float64(T), residual, contour)
end

"""
Ion-sphere radius in Bohr for mass number `A` (in atomic mass units) and mass
density `rho_m` (in g/cm^3).
"""
function ion_sphere_radius_bohr(A::Real, rho_m::Real)
    A > 0 || throw(ArgumentError("A must be positive"))
    rho_m > 0 || throw(ArgumentError("rho_m must be positive"))
    rho_kg_m3 = rho_m * GRAM_PER_CM3_TO_KG_PER_M3
    R_m = (3 * A * ATOMIC_MASS_UNIT_KG / (4π * rho_kg_m3))^(1 / 3)
    return R_m / BOHR_RADIUS_M
end

# Precompute what the DVR cardinal interpolation needs on the l = 0 grid.
function _dvr_interpolation_setup(dvr::DVR)
    N = length(dvr.grid)
    norms = [jacobi_normsq(n, 0, 0) for n = 0:(N-1)]
    return (dvr = dvr, norms = norms, transform = transpose(dvr.transformation))
end

"""
Evaluate the DVR cardinal interpolant of `values` at a single radius `r`.

The interpolant is a polynomial of degree `N - 1` through the l = 0
Gauss-Radau nodes. `values` holds a smooth radial function (the smooth part of
the effective potential), so this is spectrally accurate in the cell.
"""
function _interpolate_dvr(values::AbstractVector, setup, r::Real, R::Real)
    N = length(setup.norms)
    x = 2r / R - 1
    w = jacobi_weight(x, 0, 0)
    P = vec(jacobi_polynomials([x], 0, 0, N))
    φ = (sqrt(w) ./ sqrt.(setup.norms)) .* P
    π = setup.transform * φ
    coefficients = setup.dvr.coefficients .* values
    return sum(coefficients .* π)
end

# Antiderivative of the l = 0 FBR basis at a single point:
# `I[n] = ∫_{-1}^{x} φ_n(t) dt`, using the Legendre identity
# ∫ P_n = (P_{n+1} - P_{n-1}) / (2n + 1).
function _antiderivative_row(x::Real, N::Integer)
    P = jacobi_polynomials([x], 0, 0, N + 2)
    I = Vector{Float64}(undef, N)
    I[1] = (x + 1) / sqrt(2.0)
    for n = 1:(N-1)
        I[n+1] = (P[n+2, 1] - P[n, 1]) / (2n + 1) / sqrt(2 / (2n + 1))
    end
    return I
end

"""
Return `c(r) = ∫_0^r f(r') dr'` for values `f` sampled on the l = 0 grid.

The integration is performed on the polynomial interpolant in the FBR basis,
so for a smooth density it is accurate to roundoff.
"""
function _cumulative_integral(values::AbstractVector, dvr::DVR, R::Real)
    N = length(dvr.grid)
    coefficients = dvr.fbr_basis * (dvr.weights .* values)
    return function (r::Real)
        x = 2r / R - 1
        return (R / 2) * sum(_antiderivative_row(x, N) .* coefficients)
    end
end

"""
Hartree potential `V_H(r) = 4π[1/r ∫_0^r r'² n_e dr' + ∫_r^R r' n_e dr']`.

`density` is sampled on `x_to_r.(dvr.grid, R)`, the l = 0 common grid.
"""
function _hartree_potential(density::AbstractVector, dvr::DVR, R::Real)
    r = x_to_r.(dvr.grid, R)
    g = 4π .* (r .^ 2) .* density
    h = 4π .* r .* density
    cumulative_g = _cumulative_integral(g, dvr, R)
    cumulative_h = _cumulative_integral(h, dvr, R)
    total_h = cumulative_h(R)
    return [
        cumulative_g(r[i]) / r[i] + (total_h - cumulative_h(r[i])) for i in eachindex(r)
    ]
end

function _effective_potential(W::AbstractVector, setup, R::Real, Z::Real)
    return function (r::Real)
        r >= R && return 0.0
        return -Z / r + _interpolate_dvr(W, setup, r, R)
    end
end

"""
Self-consistent average-atom result.

Fields:
- `density`: Electron density `n_e(r)` in bohr^-3 on `grid`
- `potential`: Final effective potential `V_eff(r)` in Hartree on `grid`
- `chemical_potential_eV`: Chemical potential `μ` in eV
- `radius_bohr`: Ion-sphere radius `R` in Bohr
- `temperature_eV`: Temperature in eV
- `neutrality_residual`: `4π ∫_0^R r² n_e dr - Z`
- `self_consistency_residual`: `max_r |V_eff,next(r) - V_eff(r)|` in Hartree,
  where `V_eff,next` is the mixed update of SPEC step 5
- `iterations`: Number of self-consistent iterations performed
- `converged`: Whether both convergence criteria were met
- `states`: Normalized Siegert states for each angular momentum
- `grid`: Common radial grid in Bohr (the l = 0 Gauss-Radau grid)
- `contour`: Contour parameters used for the final density
"""
struct AverageAtomResult
    density::Vector{Float64}
    potential::Vector{Float64}
    chemical_potential_eV::Float64
    radius_bohr::Float64
    temperature_eV::Float64
    neutrality_residual::Float64
    self_consistency_residual::Float64
    iterations::Int
    converged::Bool
    states::Vector{NormalizedSPS}
    grid::Vector{Float64}
    contour::Contour
end

"""
Internal self-consistent driver with `R` in Bohr and `T` in Hartree.

See [`average_atom`](@ref) for the public interface and units.
"""
function _average_atom_bohr(
    Z::Real,
    R::Real,
    T::Real;
    N::Integer = 50,
    l_max::Integer = 10,
    xc = nothing,
    mixing::Real = 0.4,
    max_iterations::Integer = 100,
    scf_tolerance::Real = 1e-4,
    free_electron_subtraction::Bool = false,
    E_min::Real = 0.0,
    E_max::Union{Nothing,Real} = nothing,
    E_I::Union{Nothing,Real} = nothing,
    n_vert::Integer = 1,
    n_horiz::Integer = 20,
    n_return::Integer = 20,
    bound_tolerance::Real = 1e-8,
    mu_tolerance::Real = 1e-10,
)
    Z > 0 || throw(ArgumentError("Z must be positive"))
    R > 0 || throw(ArgumentError("R must be positive"))
    T > 0 || throw(ArgumentError("T must be positive"))
    N > 0 || throw(ArgumentError("N must be positive"))
    l_max >= 0 || throw(ArgumentError("l_max must be non-negative"))
    0 < mixing <= 1 || throw(ArgumentError("mixing must lie in (0, 1]"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    scf_tolerance > 0 || throw(ArgumentError("scf_tolerance must be positive"))

    dvr0 = DVR(N, 0)
    setup = _dvr_interpolation_setup(dvr0)
    r = x_to_r.(dvr0.grid, R)
    W = fill(Float64(Z) / R, length(r))
    mu_guess = nothing
    aa = nothing
    iterations = 0
    converged = false
    resid_v = Inf
    resid_n = Inf
    W_used = W

    for it = 1:max_iterations
        iterations = it
        W_used = W
        V = _effective_potential(W_used, setup, R, Z)
        aa = fixed_potential_density(
            V,
            Z,
            R,
            T;
            N = N,
            l_max = l_max,
            E_min = E_min,
            E_max = E_max,
            E_I = E_I,
            n_vert = n_vert,
            n_horiz = n_horiz,
            n_return = n_return,
            bound_tolerance = bound_tolerance,
            mu_tolerance = mu_tolerance,
            mu_guess = mu_guess,
            free_electron_subtraction = free_electron_subtraction,
        )
        mu_guess = aa.chemical_potential
        V_H = _hartree_potential(aa.density, dvr0, R)
        xc_values = xc === nothing ? zeros(length(r)) : xc(aa.density, r)
        length(xc_values) == length(r) ||
            throw(DimensionMismatch("xc must return one value per grid point"))
        xc_values = Float64.(xc_values)
        xc_at_R = xc_values[end]
        W_new = V_H .+ xc_values .- xc_at_R
        delta = maximum(abs.(W_new .- W_used))
        resid_v = mixing * delta
        resid_n = abs(aa.neutrality_residual) / Z
        if resid_v < scf_tolerance && resid_n < 1e-8
            converged = true
            break
        end
        W = (1 - mixing) .* W_used .+ mixing .* W_new
    end

    if !converged
        @warn "average_atom did not converge" Z = Z R = R T = T iterations = iterations resid_v =
            resid_v resid_n = resid_n
    end

    potential = [-Z / r[i] + W_used[i] for i in eachindex(r)]
    return AverageAtomResult(
        aa.density,
        potential,
        aa.chemical_potential * HARTREE_TO_EV,
        Float64(R),
        Float64(T) * HARTREE_TO_EV,
        aa.neutrality_residual,
        resid_v,
        iterations,
        converged,
        aa.states,
        aa.grid,
        aa.contour,
    )
end

"""
    average_atom(Z, A, rho_m, T; N=50, l_max=10, xc=nothing, mixing=0.4,
                 max_iterations=100, free_electron_subtraction=false, ...)
        -> AverageAtomResult

Self-consistent finite-temperature average atom for nuclear charge `Z`, mass
number `A` (atomic mass units), mass density `rho_m` (g/cm^3), and temperature
`T` (eV).

The effective potential is `V_eff(r) = -Z/r + V_H(r) + V_xc[n](r) - V_xc[n(R)]`
and is iterated with fixed-point mixing. The Hartree potential uses the M1
contour + Matsubara density. `xc = nothing` means zero exchange-correlation.

`xc` may be a callable `xc(n::AbstractVector, r::AbstractVector)` returning the
unsubtracted `V_xc[n](r)` on the grid; the boundary subtraction is applied by
the driver. Density is returned in bohr^-3, potential in Hartree, and
`chemical_potential_eV` in eV.

`free_electron_subtraction = true` evaluates each partial wave as
`n_l - n_l^free` and adds the analytic uniform free density back (S15
Eqs. (13)/(14)). At a fixed potential it reproduces the no-subtraction density
at a larger `l_max` to ~1e-7 relative L2. The self-consistent loop with the
subtraction is stable for small `l_max` (e.g. 4-6) but the free tail is not
accurate enough at `l_max = 10` for strongly structured potentials, so the
flag is off by default.

Convergence requires `self_consistency_residual < scf_tolerance` (default
`1e-4` Ha; the M1 density's Float64 noise floor is ~1e-5) and relative
neutrality `< 1e-8`. If convergence is not reached, a warning is emitted and
`converged = false`.

# Example
```julia
result = average_atom(5.0, 10.0, 10.0, 0.2)
result.converged
result.chemical_potential_eV
result.density
```
"""
function average_atom(
    Z::Real,
    A::Real,
    rho_m::Real,
    T::Real;
    N::Integer = 50,
    l_max::Integer = 10,
    xc = nothing,
    mixing::Real = 0.4,
    max_iterations::Integer = 100,
    free_electron_subtraction::Bool = false,
    scf_tolerance::Real = 1e-4,
    E_min::Real = 0.0,
    E_max::Union{Nothing,Real} = nothing,
    E_I::Union{Nothing,Real} = nothing,
    n_vert::Integer = 1,
    n_horiz::Integer = 20,
    n_return::Integer = 20,
    bound_tolerance::Real = 1e-8,
    mu_tolerance::Real = 1e-10,
)
    T > 0 || throw(ArgumentError("T must be positive"))
    R = ion_sphere_radius_bohr(A, rho_m)
    T_Ha = T / HARTREE_TO_EV
    return _average_atom_bohr(
        Z,
        R,
        T_Ha;
        N = N,
        l_max = l_max,
        xc = xc,
        mixing = mixing,
        max_iterations = max_iterations,
        scf_tolerance = scf_tolerance,
        free_electron_subtraction = free_electron_subtraction,
        E_min = E_min,
        E_max = E_max,
        E_I = E_I,
        n_vert = n_vert,
        n_horiz = n_horiz,
        n_return = n_return,
        bound_tolerance = bound_tolerance,
        mu_tolerance = mu_tolerance,
    )
end
