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
# Below are two simple illustrations.
#
#
#
# ### Finite region with non-zero boundary
#
# First of all, they are often defined only in some region ``Ω``, and do not zero out on the boundary of the region ``∂Ω```.
# This is simulated below by adding some linear tilt to the phase.

tilt = [(i[1] * 0.000135 + i[2] * 0.0011 + 10) for i in eachindex(IndexCartesian(), mask)]
resp_wr = phwrap((resp .+ tilt)) .* mask
phasedisplay(resp_wr .* ap2mask(mask))
# (here we have used [`ap2mask`](@ref) function which removes the pixels outside the mask for better visibility).

# Here is what happens if we subtract the (exactly known) tilt from our new wrapped phase

phasedisplay((resp_wr .- tilt) .* ap2mask(mask))
phasedisplay(resp_wr .* ap2mask(mask))
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

# ### Small noise
#
# The straightforward unwrapping becomes impossible if value of one pixel has jumped to approximately on π radians (for instance, due to a broken pixel).
# You can hardly see the difference (check pixel (100,100)).
addpi!(arr, r, c) = (arr[r, c] = phwrap(arr[r, c] + π))

resp_wr = phwrap(resp)
addpi!(resp_wr, 100, 100)
phasedisplay(resp_wr .* ap2mask(mask))

# Try to unwrap it:
ph_itoh = reshape(itoh(resp_wr[:]), size(resp_wr))
arraydisplay((ph_itoh) .* mask)

# and we see that one pixel made the unwrapping impossible, because a spurious phase jump on 2π appeared on passing the broken pixel.
# But we can unwrap it with the least-squares algorithm:
resp_LS = unwrap_LS(resp_wr, mask; restore_piston=true)
fig, ax, hm = arraydisplay(resp_LS);
Colorbar(fig[1, 2], hm)
fig













# ## Example with wedge (with zero mean for simplicity of comparison with the restored)
wedge = mask .* linearphase((300, 300), 150, 150, 0.2, 0.11)
piston = maskedrmse(wedge, mask)
wedge .-= piston
fig, ax, hm = arraydisplay(wedge)
Colorbar(fig[1, 2], hm)
ax.title = "Original wedge"
fig

wedge_wr = phwrap(wedge)
phasedisplay(wedge_wr)

# It's impossible to unwrap it with the Itoh algorithm.
wedge_itoh = reshape(itoh(wedge_wr[:]), size(wedge_wr))
fig, ax, hm = arraydisplay(wedge_itoh)
Colorbar(fig[1, 2], hm)
ax.title = "RMS error is $(maskedrmse(wedge_itoh, wedge, mask))"
fig

# Nor it's possible to do this using Itoh's algorithm along the rows.
wedge_itoh = reshape(itoh(wedge_wr'[:]), size(wedge_wr))'
fig, ax, hm = arraydisplay(wedge_itoh)
Colorbar(fig[1, 2], hm)
ax.title = "RMS error is $(maskedrmse(wedge_itoh, wedge, mask))"
fig

# But we can unwrap it with the least-squares algorithm (up to an integer multiple of 2π, which we calculate here as the value of the error at the central pixel)
wedge_LS = unwrap_LS(wedge_wr, mask; restore_piston=false)
piston = wedge_LS[150, 150] - wedge[150, 150]
wedge_LS .-= piston
fig, ax, hm = arraydisplay(wedge_LS);
Colorbar(fig[1, 2], hm)
fig

# Compare with the original wedge:
fig, ax, hm = arraydisplay((wedge_LS .- wedge))
Colorbar(fig[1, 2], hm)
ax.title = "RMS error is $(maskedrmse(wedge_LS, wedge, mask))"
fig

# ## Unwrapping of the phase with residues
# Let us add some bad pixel groups to the wrapped wedge
resblock!(arr, r, c) = (arr[c:(c + 1), r:(r + 1)] .= [0 π/2; -π/2 π])
resblock!(wedge_wr, 100, 110)
resblock!(wedge_wr, 50, 150)
resblock!(wedge_wr, 130, 130)

# And add 1 π to some pixels
addpi!(wedge_wr, 200, 200)
addpi!(wedge_wr, 40, 180)

fig, ax, hm = phasedisplay(wedge_wr)
ax.title = "RMS error with the original phase is $(maskedrmse(wedge_wr, phwrap(wedge), mask))"
fig

# It's still impossible to unwrap it with the Itoh algorithm, but now the error is of a different type
wedge_itoh = reshape(itoh(wedge_wr[:]), size(wedge_wr))
fig, ax, hm = arraydisplay(wedge_itoh)
ax.title = "RMS error is $(maskedrmse(wedge_itoh, wedge, mask))"
Colorbar(fig[1, 2], hm)
fig


# Nor it's possible to do this using Itoh's algorithm along the rows.
wedge_itoh = reshape(itoh(wedge_wr'[:]), size(wedge_wr))'
fig, ax, hm = arraydisplay(wedge_itoh)
ax.title = "RMS error is $(maskedrmse(wedge_itoh, wedge, mask))"
Colorbar(fig[1, 2], hm)
fig

# But we can unwrap it with the least-squares algorithm
wedge_LS = unwrap_LS(wedge_wr, mask; restore_piston=true)
wedge_LS .-= piston
fig, ax, hm = arraydisplay(wedge_LS)
Colorbar(fig[1, 2], hm)
fig

# Compare with the original wedge:
fig, ax, hm = arraydisplay(ap2mask(mask) .* (wedge_LS .- wedge))
ax.title = "RMS error is $(maskedrmse(wedge_LS, wedge, mask))"
Colorbar(fig[1, 2], hm)
fig
