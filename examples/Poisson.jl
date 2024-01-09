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

mask = zeros(300, 300)
mask = [
    ((i[1] - 150)^2 + (i[2] - 150)^2 < 145^2) for i in eachindex(IndexCartesian(), mask)
]
act = [((i[1] - 50)^2 + (i[2] - 153.5)^2 < 10^2) for i in eachindex(IndexCartesian(), mask)];

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

# # Phase Unwrapping
# We will use the calculated response function for the phase unwrapping demonstration.
# First we define a couple of plotting functions to save on typing
phasedisplay(args...; kwargs...) = heatmap(
    args...;
    axis=(aspect=DataAspect(),),
    colormap=:cyclic_mygbm_30_95_c78_n256,
    kwargs...,
)
arraydisplay(args...; kwargs...) = heatmap(args...; axis=(aspect=DataAspect(),), kwargs...)

resp ./= 10
arraydisplay(resp)

resp_wr = phwrap(resp)
phasedisplay(resp_wr)

# This noiseless, periodic, and defined on the whole domain phase is easy to unwrap by line-by-line integration of wrapped phase differences using Itoh algorithm, see [`itoh`](@ref).
ph_itoh = reshape(itoh(resp_wr[:]), size(resp_wr))
arraydisplay(ph_itoh)

# Phases occurring in real life are usually a little bit more complicated.
# First of all, they are often defined only in some region ``Ω```, and do not zero out on the boundary of the region ``∂Ω```.
# This is simulated below by adding some linear tilt to the phase.

tilt = [(i[1] * 0.000135 + i[2] * 0.0011 + 10) for i in eachindex(IndexCartesian(), mask)]
resp_wr = phwrap((resp .+ tilt)) .* mask
phasedisplay(resp_wr .* ap2mask(mask))
# (here we have used [`ap2mask`](@ref) function which removes the pixels outside the mask for better visibility).

# Here is what happens if we subtract the (exactly known) tilt from our new wrapped phase

phasedisplay((resp_wr .- tilt) .* ap2mask(mask))

#

#
# This phase should not be restorable by the Itoh algorithm because the jumps on the boundary are greater than ``π``, and indeed we see the errors in unwrapping
ph_itoh = reshape(itoh(resp_wr[:]), size(resp_wr))
arraydisplay((ph_itoh .- tilt) .* mask)
# but wrapped back it looks OK, of course (what we see is just unwrapping errors, assigning wrong ``k`` to an unwrapped value `` \hat{\phi} = \psi + 2πk, k ∈ \mathbb{Z}``):
phasedisplay(phwrap(ph_itoh .- tilt) .* ap2mask(mask))

# We can try to unwrap it with the least-squares algorithm
phi_LS = unwrap_LS(resp_wr, mask)
arraydisplay((phi_LS .- tilt) .* mask)
# The error here is more fundamental:
phasedisplay(phwrap(phi_LS .- tilt) .* ap2mask(mask))
