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

# ## Conclusions
#
# It works as expected.
#
#
#
#
#
#
#
#
# # Dual and half-dual grids
#
# For Talmi-Ribak's method, I also need to find points on the dual grid that fall within the area and on its boundary (with at least two adjacent pixels).
# Here is a straightforward approach.
#
function dual_region(ap)
    dualap = zeros(Bool, size(ap) .+ 1) # index i,j corresponds to the point i-1/2,j-1/2
    for ind in eachindex(IndexCartesian(), ap)
        if ap[ind]
            if ap[ind + CartesianIndex(1, 0)] #two horisontal neighbours (i,j) and (i+1,j)
                dualap[ind + CartesianIndex(1, 0)] = 1 # (i+1/2,j-1/2)
                dualap[ind + CartesianIndex(1, 1)] = 1 # (i+1/2,j+1/2)
            end
            if ap[ind + CartesianIndex(0, 1)] #two vertical neighbours (i,j) and (i,j+1)
                dualap[ind + CartesianIndex(0, 1)] = 1 # (i-1/2,j+1/2)
                dualap[ind + CartesianIndex(1, 1)] = 1 # (i+1/2,j+1/2)
            end
        end
    end
    return dualap
end

function dualcoordmin(ind)
    return ind.I .- 0.5
end
function dualcoordplus(ind)
    return ind.I .+ 0.5
end
#
# Plot the ap and its dual region.
#
fig, ax, hm = arraydisplay(ap; colormap=:reds, alpha=0.25);
dualap = dual_region(ap)
cwborder, _ = PhaseUtils._find_cw_border_alloc(dualap)
scatter!(dualcoordmin.(findall(dualap)); marker='+', color=:blue, label="dual region")
scatter!(map(dualcoordmin, cwborder); color=:blue, marker='□', label="dual contour")
leg = axislegend(ax)
ax.title = "Detected dual region and its border"
fig

# # Alternative approach
#
# Alternative approach would be to look where the residue can be defined, and this mean all 4 pixels surroundg the dual vertex should be in the aperture.
# This is equivalent to covolution  with a  box 2x2.

function dual_region_box(ap)
    dualap = zeros(Bool, size(ap) .+ 1) # index i,j corresponds to the point i-1/2,j-1/2
    # dualap[1:(end - 1), 1:(end - 1)] .= ap
    for ind in eachindex(IndexCartesian(), dualap[2:(end - 1), 2:(end - 1)])
        dualap[ind] =
            ap[ind] &&
            ap[ind - CartesianIndex(1, 1)] &&
            ap[ind - CartesianIndex(1, 0)] &&
            ap[ind - CartesianIndex(0, 1)]
    end
    return dualap
end

# Plot the ap and its dual region.
#
fig, ax, hm = arraydisplay(ap; colormap=:reds, alpha=0.25);
dualap = dual_region_box(ap)
cwborder, _ = PhaseUtils._find_cw_border_alloc(dualap; outside=true)
scatter!(dualcoordmin.(findall(dualap)); marker='+', color=:blue, label="dual region")
scatter!(map(dualcoordmin, cwborder); color=:blue, marker='□', label="dual contour")
leg = axislegend(ax)
ax.title = "Detected dual region and its border"
fig
