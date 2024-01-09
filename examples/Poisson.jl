# ```@meta
# CurrentModule = PhaseUtils
# DocTestSetup = quote
#     using PhaseUtils
# end
# ```

using PhaseUtils
using CairoMakie
CairoMakie.activate!(; type="png")

# # Solving Poisson equation

# Let us simulate a response of a membrane deformable mirror using the Poisson equation.
# Define a circular aperture of 145-pixel radius and a circular actuator  of 10-pixel radius inside it:

mask = zeros(300, 300)
mask = [
    ((i[1] - 150)^2 + (i[2] - 150)^2 < 145^2) for i in eachindex(IndexCartesian(), mask)
]
act = [((i[1] - 120)^2 + (i[2] - 150)^2 < 10^2) for i in eachindex(IndexCartesian(), mask)]

# Plot actuator inside the aperture:
heatmap(mask .+ act; axis=(aspect=DataAspect(),))

# Calculate the response and show the results
resp = membrane_sor(act, mask)
contourf(resp; axis=(aspect=DataAspect(),))

# Check that the Laplacian of the calculated response is proportional to the actuator shape:
act_restored = PhaseUtils._calculate_Laplacian(resp)
heatmap(act_restored .* mask; axis=(aspect=DataAspect(),))


# To set the boundary conditions, use the mutating version of the solver:
mask = zeros(Bool, 300, 300)
mask[2:(end - 1), 2:(end - 1)] .= true
u = zeros(size(mask))
u[1, :] .= 1
membrane_sor!(u, zeros(size(u)), mask)
heatmap(u; axis=(aspect=DataAspect(),))
