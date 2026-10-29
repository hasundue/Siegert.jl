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
Average-atom Siegert state.

Wraps a raw `SPS` and adds the data needed by the average-atom layer: the DVR
basis, the normalized coefficients, and the raw normalization constants.

Fields:
- `sps`: Raw eigenproblem result
- `coefficients`: Normalized DVR coefficients
- `normalizations`: Raw normalization constants `N_n`
- `dvr`: DVR basis used to build and evaluate the state
"""
struct AverageAtomSPS <: AbstractSPS
    sps::SPS
    coefficients::Matrix{ComplexF64}
    normalizations::Vector{ComplexF64}
    dvr::DVR
end

# Accessors shared with the raw SPS interface.
momenta(aa::AverageAtomSPS) = momenta(aa.sps)
coefficients(aa::AverageAtomSPS) = aa.coefficients
grid(aa::AverageAtomSPS) = grid(aa.sps)
box_size(aa::AverageAtomSPS) = box_size(aa.sps)
angular_momentum(aa::AverageAtomSPS) = angular_momentum(aa.sps)
hamiltonian(aa::AverageAtomSPS) = hamiltonian(aa.sps)
metric(aa::AverageAtomSPS) = metric(aa.sps)
boundary(aa::AverageAtomSPS) = boundary(aa.sps)
dvr(aa::AverageAtomSPS) = aa.dvr

"""
Compute the SPS normalization constants and rescale the DVR coefficients.

The normalization follows SS23 Eq. (6) and BT07 Eq. (37). The radial
wavefunction is `P_n(r) = r * Σ_i c_{i,n} π_i(x)`, so the volume term is
`(a/2) c_n^T ρ c_n` and the surface value is `P_n(a) = a π(a)^T c_n`.

Returns the raw normalization constants `N_n`.
"""
function _normalize_average_atom!(
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
    AverageAtomSPS(sps::SPS) -> AverageAtomSPS
    AverageAtomSPS(sps::SPS, dvr::DVR) -> AverageAtomSPS

Wrap a raw `SPS` and apply SS23 Eq. (6) / BT07 Eq. (37) normalization. The
optional `dvr` avoids rebuilding the basis when it is already available.
"""
function AverageAtomSPS(sps::SPS, basis::DVR)
    N = size(sps.coefficients, 1)
    size(basis.basis, 1) == N ||
        throw(DimensionMismatch("DVR basis does not match the SPS dimension"))
    coefficients = copy(sps.coefficients)
    norm = _normalize_average_atom!(
        coefficients,
        sps.momenta,
        sps.metric,
        basis,
        sps.box_size,
        sps.angular_momentum,
    )
    return AverageAtomSPS(sps, coefficients, norm, basis)
end

function AverageAtomSPS(sps::SPS)
    basis = DVR(size(sps.coefficients, 1), sps.angular_momentum)
    return AverageAtomSPS(sps, basis)
end

"""
Result of a fixed-potential average-atom density calculation.

Fields:
- `density`: Electron density `n_e(r)` on `grid`
- `chemical_potential`: Chemical potential `μ` in Hartree
- `states`: Siegert states for each angular momentum
- `grid`: Common radial grid
- `radius`: Ion-sphere radius `R`
- `temperature`: Temperature in Hartree
- `neutrality_residual`: `4π ∫ r² n_e dr - Z` at the returned `μ`
- `contour`: Contour parameters used for the final density
"""
struct AADensity
    density::Vector{Float64}
    chemical_potential::Float64
    states::Vector{AverageAtomSPS}
    grid::Vector{Float64}
    radius::Float64
    temperature::Float64
    neutrality_residual::Float64
    contour::Contour
end

"""
Return the `l = 0` Gauss-Radau grid and the weights for `∫_0^R f(r) dr`.

The same grid is used for every angular momentum because the DVR grids depend
on `l`. For an `AverageAtomSPS` with `l = 0` its own grid is the common grid.
"""
function common_radial_grid(aa::AverageAtomSPS)
    N = size(coefficients(aa), 1)
    basis = angular_momentum(aa) == 0 ? dvr(aa) : DVR(N, 0)
    r = x_to_r.(basis.grid, box_size(aa))
    weights = (box_size(aa) / 2) .* basis.weights
    return r, weights
end

function common_radial_grid(states::AbstractVector{<:AverageAtomSPS})
    isempty(states) && throw(ArgumentError("states must not be empty"))
    return common_radial_grid(first(states))
end

"""
Return the weights for the volume integral `∫_0^R r² f(r) dr`.
"""
volume_weights(r::AbstractVector{<:Real}, weights::AbstractVector{<:Real}) =
    (r .^ 2) .* weights

"""
    siegert_wavefunction(aa::AverageAtomSPS) -> Matrix{ComplexF64}

Evaluate all normalized radial Siegert states of `aa` on `grid(aa)`.

The returned matrix has one row per grid point and one column per state, with
`P_{n,l}(r) = r Σ_i c_{i,n} π_i(x(r))`.
"""
function siegert_wavefunction(aa::AverageAtomSPS)
    return siegert_wavefunction(aa, grid(aa))
end

siegert_wavefunction(sps::SPS) =
    throw(ArgumentError("wrap the raw SPS in AverageAtomSPS first"))
siegert_wavefunction(sps::SPS, r::AbstractVector) =
    throw(ArgumentError("wrap the raw SPS in AverageAtomSPS first"))

"""
    siegert_wavefunction(aa::AverageAtomSPS, r::AbstractVector) -> Matrix{ComplexF64}

Evaluate all normalized radial Siegert states of `aa` on the radii `r`.

`r` must lie in `[0, R]`. The returned matrix has one row per radius and one
column per state.
"""
function siegert_wavefunction(aa::AverageAtomSPS, r::AbstractVector)
    all(0 .<= r .<= box_size(aa)) || throw(DomainError(r, "radii must lie in [0, R]"))
    l = angular_momentum(aa)
    C = coefficients(aa)
    N = size(C, 1)
    x = 2 .* r ./ box_size(aa) .- 1

    w = jacobi_weight.(x, 0, 2l)
    h = jacobi_normsq.(0:(N-1), 0, 2l)
    P = jacobi_polynomials(x, 0, 2l, N)
    φ = sqrt.(reshape(w, 1, :) ./ reshape(h, :, 1)) .* P
    π = transpose(dvr(aa).transformation) * φ
    values = transpose(C) * π
    return ComplexF64.(transpose(reshape(r, 1, :) .* values))
end

"""
    siegert_green(aa::AverageAtomSPS, k::Number; r=grid(aa)) -> Vector{ComplexF64}

Diagonal partial-wave Green's function
`G_l(r,r;k) = Σ_n P_{n,l}(r)^2 / [k_n (k - k_n)]` on the radii `r`.
"""
function siegert_green(aa::AverageAtomSPS, k::Number; r = grid(aa))
    P = siegert_wavefunction(aa, r)
    return _green_diagonal(P, momenta(aa), k)
end

function _green_diagonal(P::AbstractMatrix, ks::AbstractVector, k::Number)
    state_weights = @. 1 / (ks * (k - ks))
    return (P .^ 2) * state_weights
end

"""
    siegert_dos(aa::AverageAtomSPS, k::Number; r, weights) -> Float64

Partial density of states `χ_l(k) = -(1/π) Im ∫_0^R dr 2(2l+1) G_l(r,r;k)`.

When `r` and `weights` are omitted, the common `l = 0` grid is used. The
weights must integrate over `dr`, not over `r² dr`.
"""
function siegert_dos(
    aa::AverageAtomSPS,
    k::Number;
    r::Union{Nothing,AbstractVector} = nothing,
    weights::Union{Nothing,AbstractVector} = nothing,
)
    if r === nothing
        r, default_weights = common_radial_grid(aa)
        weights === nothing && (weights = default_weights)
    end
    weights === nothing && throw(ArgumentError("weights must be given together with r"))
    length(r) == length(weights) ||
        throw(DimensionMismatch("r and weights must have the same length"))
    G = siegert_green(aa, k; r = r)
    return real(-(1 / π) * imag(2 * (2angular_momentum(aa) + 1) * sum(weights .* G)))
end

"""
    siegert_dos(states::AbstractVector{<:AverageAtomSPS}, k::Number; r, weights)
    -> Float64

Total density of states, summed over the given partial waves.
"""
function siegert_dos(
    states::AbstractVector{<:AverageAtomSPS},
    k::Number;
    r::Union{Nothing,AbstractVector} = nothing,
    weights::Union{Nothing,AbstractVector} = nothing,
)
    if r === nothing
        r, default_weights = common_radial_grid(states)
        weights === nothing && (weights = default_weights)
    end
    weights === nothing && throw(ArgumentError("weights must be given together with r"))
    return sum(siegert_dos(aa, k; r = r, weights = weights) for aa in states)
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

function _state_data(states::AbstractVector{<:AverageAtomSPS}, r::AbstractVector)
    return [
        (momenta(aa), siegert_wavefunction(aa, r) .^ 2, angular_momentum(aa)) for
        aa in states
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
    return density, contour
end

"""
    siegert_density(states, μ, T; E_min=0, r, weights, ...) -> Vector{Float64}

Electron density `n_e(r)` from SS23 Eq. (13) on the common radial grid.

Bound states with energy below `E_min` are added explicitly; the continuum is
evaluated with the three-sided contour plus the enclosed Matsubara poles.
"""
function siegert_density(
    states::AbstractVector{<:AverageAtomSPS},
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
    )
    return density
end

"""
    chemical_potential(states, Z, R, T; E_min=0, ...) -> (μ, density, residual)

Solve `Z = 4π ∫_0^R r² n_e(r) dr` for the chemical potential by bisection with
an expanding bracket. The residual is `4π ∫ r² n_e dr - Z`.
"""
function chemical_potential(
    states::AbstractVector{<:AverageAtomSPS},
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
        )
        electron_number = 4π * sum(volume .* density)
        return electron_number - Z, density, contour
    end

    energies = reduce(vcat, (real.(momenta(aa) .^ 2 ./ 2) for aa in states))
    lower = minimum(energies) - 10T
    upper = maximum(energies) + 10T
    f_lower, _, _ = evaluate(lower)
    f_upper, _, _ = evaluate(upper)

    scale = max(1.0, upper - lower)
    iterations = 0
    while f_lower > 0 && iterations < max_iterations
        lower -= scale
        scale *= 2
        f_lower, _, _ = evaluate(lower)
        iterations += 1
    end
    scale = max(1.0, upper - lower)
    while f_upper < 0 && iterations < max_iterations
        upper += scale
        scale *= 2
        f_upper, _, _ = evaluate(upper)
        iterations += 1
    end
    (f_lower <= 0 && f_upper >= 0) ||
        throw(ArgumentError("could not bracket the chemical potential"))

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
)
    N > 0 || throw(ArgumentError("N must be positive"))
    l_max >= 0 || throw(ArgumentError("l_max must be non-negative"))
    states = map(0:l_max) do l
        sps, basis = _solve_sps(N, l, V, R)
        AverageAtomSPS(sps, basis)
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
    )
    grid, _ = common_radial_grid(states)
    contour = contour_parameters(mu, T; E_min = E_min, E_max = E_max, E_I = E_I)
    return AADensity(density, mu, states, grid, Float64(R), Float64(T), residual, contour)
end
