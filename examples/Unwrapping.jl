# ```@meta
# CurrentModule = PhaseUtils
# DocTestSetup = quote
#     using PhaseUtils
# end
# ```
#

# # Phase Unwrapping
# We will use the calculated response function for the phase unwrapping demonstration.

using PhaseUtils
using CairoMakie
CairoMakie.activate!(; type="png")
mask = circlemask((300, 300), 150, 150, 145)
act = circlemask((300, 300), 50, 153.5, 10)
resp = membrane_sor(act, mask);

# First, we define a couple of plotting functions to save on typing
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

# This noiseless, periodic, and defined on the whole domain phase is easy to unwrap by line-by-line integration of wrapped phase differences using the Itoh algorithm, see [`itoh`](@ref).
ph_itoh = reshape(itoh(resp_wr[:]), size(resp_wr))
arraydisplay(ph_itoh)

# Phases occurring in real life are usually a little bit more complicated.
# First of all, they are often defined only in some region ``Ω``, and do not zero out on the boundary of the region ``∂Ω```.
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

# We can unwrap it with the least-squares algorithm
phi_LS = unwrap_LS(resp_wr, mask)
arraydisplay((phi_LS .- tilt) .* ap2mask(mask))
