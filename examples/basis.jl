using Plots
using Siegert

# Example parameters
N = 8
l = 1
α = 0.0
β = 2l

# 1) Orthonormal L2 Jacobi basis (α, β)
φ = jacobi_basis_L2(N, α, β)
plt1 = plot(
    xlabel = "x",
    ylabel = "φ_n(x) = √w(x) P_n/√h_n",
    legend = :topright,
    title = "Jacobi L2 basis (α=$α, β=$β)",
)
for n = 1:N
    plot!(plt1, φ[n], -1, 1, label = "n=$n")
end
plot!(plt1, xlims = (-1, 1))
display(plt1)

# 2) Jacobi Gaussian DVR basis (l-based API)
ψ, z, Λ = jacobi_dvr_basis(N, l)
plt2 = plot(
    xlabel = "x",
    ylabel = "ψ_i(x)",
    legend = :topright,
    title = "Jacobi DVR basis (l=$l)",
)
for i = 1:N
    plot!(plt2, ψ[i], -1, 1, label = "i=$i")
end
# Overlay grid points and expected cardinality scale 1/√Λ
scatter!(plt2, z, zeros(length(z)), markershape = :vline, label = "grid z_i")
scatter!(plt2, z, 1.0 ./ sqrt.(Λ), label = "1/√Λ_i", ms = 3)
plot!(plt2, xlims = (-1, 1))
display(plt2)
