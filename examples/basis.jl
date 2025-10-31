using Plots
include(joinpath(@__DIR__, "..", "src", "jacobi.jl"))

# Example parameters
N = 8
l = 0
α = 0.0
β = 2l

φ = jacobi_basis_L2(N, α, β)

plt = plot(xlabel = "x", ylabel = "φ_n(x) = √w(x) P_n/√h_n", legend = :topright)
for n = 1:N
    plot!(plt, φ[n], -1, 1, label = "n=$n")
end
plot!(plt, xlims = (-1, 1))
display(plt)
