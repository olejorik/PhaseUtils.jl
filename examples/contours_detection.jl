# ```@meta
# CurrentModule = PhaseUtils
# DocTestSetup = quote
#     using PhaseUtils
# end
# ```
#

# # Test of the contour detection
#
# I use simple contour detection for phase unwrapping.
# The goal is to find the points where the finite difference is not defined.
#
# Begin with a circular aperture.
using PhaseUtils
using CairoMakie
CairoMakie.activate!(; type="png")
arraydisplay(args...; kwargs...) = heatmap(args...; axis=(aspect=DataAspect(),), kwargs...)

ap = circlemask((30, 20), 10, 10, 8)
fig, ax, hm = arraydisplay(ap; colormap=:reds, alpha=0.25)
edg = PhaseUtils._find_set_edges(ap)
for k in keys(edg)
    scatter!(
        map(x -> Tuple(x) .+ Tuple(PhaseUtils.step[k]) .* 0.2, edg[k]); label="$(String(k))"
    )
end
leg = axislegend(ax)
ax.title = "Detected edge pixels for inner circle"
fig

# We can now sort all these pixels in a clockwise manner
cwborder, _ = PhaseUtils._find_cw_border_alloc(ap)
scatter!(map(Tuple, cwborder); marker='o', color=:black, label="Contour")
delete!(leg)
axislegend(ax)
fig

# Now check the same for the inverse mask
notap = .!ap
fig, ax, hm = arraydisplay(notap; colormap=:reds, alpha=0.25)
edg = PhaseUtils._find_set_edges(notap)
for k in keys(edg)
    scatter!(
        map(x -> Tuple(x) .+ Tuple(PhaseUtils.step[k]) .* 0.2, edg[k]); label="$(String(k))"
    )
end
leg = axislegend(ax)
ax.title = "Detected edge pixels for outer circle"
fig

# We can now sort all these pixels in a clockwise manner
cwborder, _ = PhaseUtils._find_cw_border_alloc(ap; outside=true)
scatter!(map(Tuple, cwborder); marker='o', color=:black, label="Contour")
delete!(leg)
axislegend(ax)
fig

# # More difficult shapes
ap = circlemask((40, 22), 10, 11, 8) .|| circlemask((40, 22), 26, 11, 9)
fig, ax, hm = arraydisplay(ap; colormap=:reds, alpha=0.25)
edg = PhaseUtils._find_set_edges(ap)
for k in keys(edg)
    scatter!(
        map(x -> Tuple(x) .+ Tuple(PhaseUtils.step[k]) .* 0.2, edg[k]); label="$(String(k))"
    )
end
leg = axislegend(ax)
ax.title = "Detected edge pixels for inner edge"
fig

# We can now sort all this pixels in a clockwise manner
cwborder, _ = PhaseUtils._find_cw_border_alloc(ap)
scatter!(map(Tuple, cwborder); marker='o', color=:black, label="Contour")
delete!(leg)
axislegend(ax)
fig

# Now check the same for the inverse mask
notap = .!ap
fig, ax, hm = arraydisplay(notap; colormap=:reds, alpha=0.25)
edg = PhaseUtils._find_set_edges(notap)
for k in keys(edg)
    scatter!(
        map(x -> Tuple(x) .+ Tuple(PhaseUtils.step[k]) .* 0.2, edg[k]); label="$(String(k))"
    )
end
leg = axislegend(ax)
ax.title = "Detected edge pixels for outer edge"
fig

# We can now sort all these pixels in a clockwise manner
cwborder, _ = PhaseUtils._find_cw_border_alloc(ap; outside=true)
scatter!(map(Tuple, cwborder); marker='o', color=:black, label="Contour")
delete!(leg)
axislegend(ax)
fig

# # Conclusions
#
# It works as expected.
