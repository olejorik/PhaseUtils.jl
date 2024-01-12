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
#
# Use [`membrane_sor(f, a)`](@ref) to calculate the solution of the Poisson equation with zero boundary conditions
# ```math
# \begin{aligned}
# \Delta u(x) = & f(x), \quad x \in Ω,\\
# u(x) =  & 0, \quad  x \in \partial Ω.
# \end{aligned}
# ```
# Region ``Ω`` is defined by the boolean array `a`, which should have the same dimensions as the source array `f`.
#
#


# Let us simulate a response of a membrane deformable mirror using the Poisson equation.
# Define a circular aperture of 145-pixel radius and a circular actuator  of 10-pixel radius inside it:

mask = circlemask((300, 300), 150, 150, 145)
act = circlemask((300, 300), 50, 153.5, 10);

# Plot actuator inside the aperture:
heatmap(mask .+ act; axis=(aspect=DataAspect(),))

# Calculate the response and show the results
resp = membrane_sor(act, mask)
contourf(resp; axis=(aspect=DataAspect(),))

# Check that the Laplacian of the calculated response is proportional to the actuator shape:
act_restored = PhaseUtils._calculate_Laplacian(resp)
heatmap(act_restored .* mask; axis=(aspect=DataAspect(),))


# To set the non-zero boundary conditions, use the mutating version of the solver:
maskrect = zeros(Bool, 300, 300)
maskrect[2:(end - 1), 2:(end - 1)] .= true
u = zeros(size(mask))
u[1, :] .= 1
u[end, 100:200] .= -1
membrane_sor!(u, zeros(size(u)), maskrect)
fig, ax, hm = heatmap(u; axis=(aspect=DataAspect(),))
contour!(u; labels=true, levels=-1:0.1:1, labelsize=15, color=:black)
fig

# Next example sets the boundary values for a circular aperture

edge_out, _ = find_cw_border(mask; outside=true)
u = zeros(size(mask))
u[edge_out] .= sin.(range(0, 4 * 2π, length(edge_out) + 1)[1:(end - 1)])
heatmap(u; axis=(aspect=DataAspect(),))

membrane_sor!(u, zeros(size(u)), mask; maxits=1000)
fig, ax, hm = heatmap(u .* ap2mask(mask); axis=(aspect=DataAspect(),))
contour!(u; labels=true, levels=-0.9:0.1:0.9, labelsize=15, color=:white)
fig
