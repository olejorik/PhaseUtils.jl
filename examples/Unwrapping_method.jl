# ```@meta
# CurrentModule = PhaseUtils
# DocTestSetup = quote
#     using PhaseUtils
# end
# ```
#

# # Phase Unwrapping method  explained
#
# The method used in this package is based on the Least-Squares integration of the (maybe incosistent) gradient fields defiend in some connected region.

using PhaseUtils
using CairoMakie
CairoMakie.activate!(; type="png")


# # First we define several helper functions.
hmap(args...; kwargs...) =
    heatmap(args...; kwargs..., axis=merge((aspect=DataAspect(),), get(kwargs, :axis, (;))))

function residuepixel(center=[0, 0], s=(400, 400))
    x, y = 1:s[1], 1:s[2]
    return [atan(yp - center[2], xp - center[1]) for xp in x, yp in y]
end

function resphasepixel(ap, ppos, pneg)
    s = size(ap)
    ph = zeros(s)
    for p in ppos
        ph .+= residuepixel(p, s)
    end
    for p in pneg
        ph .-= residuepixel(p, s)
    end
    return phwrap(ph)
end

function show_grad(gr, label)
    fig = Figure()
    ax = CairoMakie.Axis(fig[1, 1]; aspect=DataAspect())
    hm = heatmap!(ax, gr[1])
    ax.title = "First component of the gradient $label"
    Colorbar(fig[1, 2], hm)
    ax = CairoMakie.Axis(fig[1, 3]; aspect=DataAspect())
    hm = heatmap!(ax, gr[2])
    Colorbar(fig[1, 4], hm)
    ax.title = "Second component of the gradient $label"
    rowsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)
    return fig |> display
end

function show_res!(ax, posx, posy, resmap)
    respos = resmap .> 0
    resneg = resmap .< 0
    scatter!(ax, posx[respos], posy[respos]; color="white", markersize=10)
    return scatter!(ax, posx[resneg], posy[resneg]; color="black", markersize=10)
end


function dual_region_box(ap)
    dualap = zeros(Bool, size(ap) .+ 1) # index i,j corresponds to the point i-1/2,j-1/2
    ## dualap[1:(end - 1), 1:(end - 1)] .= ap
    for ind in eachindex(IndexCartesian(), dualap[2:(end - 1), 2:(end - 1)])
        dualap[ind] =
            ap[ind] &&
            ap[ind - CartesianIndex(1, 1)] &&
            ap[ind - CartesianIndex(1, 0)] &&
            ap[ind - CartesianIndex(0, 1)]
    end
    return dualap
end

function dualcoordmin(ind)
    return ind.I .- 0.5
end
function dualcoordplus(ind)
    return ind.I .+ 0.5
end
function dualcoordplus(x, y)
    return CartesianIndex(Int(x + 0.5), Int(y + 0.5))
end

##
# We will do everything on a small resolution example phase.

n = 8
ap = circlemask((2n, 2n + 1), n, n + 1, n - 2)
mask = ap2mask(ap)
hmap(ap)

# Choose `FiniteDifferences` as default method for gradients etc
diffmeth = PhaseUtils.FiniteDifferences()
get_grad(arr) = PhaseUtils._calculate_gradient(arr, diffmeth)

# Create a test phase with a residue
phgt = linearphase(size(ap), 0, 0, 0.2, 0.7)
ph = phgt + resphasepixel(ap, [(n + 0.5, n + 1.5)], [])
phm = ph .* mask
hmap(ph; axis=(title="Original phase, not masked",)) |> display
hmap(phm; axis=(title="Origimal phase, masked",)) |> display

# Calculate gradient of the masked phase
gr = get_grad(phm)
show_grad(gr, "\nof the original phase")

# and wrapped gradient
g = (g1, g2) = phwrap.(gr)
show_grad(g, "\nwrapped")

#  Find the residues in the phase
posx, posy, resmap = getresmapsparce(phm)

show_res!(posx, posy, resmap) = show_res!(current_axis(), posx, posy, resmap)
fig, ax, hm = hmap(phm)
show_res!(posx, posy, resmap)
ax.title = "$(length(resmap)) residue(s) found"
fig |> display

# ### Build and solve the Poisson equation for the correction potential
dualap = dual_region_box(ap)
dualap_inside = copy(dualap)



rho = zeros(size(dualap_inside))
for ind in eachindex(resmap)
    rho[dualcoordplus(posx[ind], posy[ind])] = π / 2 * resmap[ind]
end
fig, ax, hm = hmap(rho);
heatmap!(ax, ap2mask(1 .- dualap_inside); colormap=:blues)
ax.title = "region and residues on the dual grid"
fig |> display

corr_pot = membrane_sor(rho, dualap_inside)
fig, ax, hm = hmap(corr_pot .* ap2mask(dualap));
ax.title = "Correcting potential"
Colorbar(fig[1, 2], hm)
fig |> display

# Calculate the gradient of the correction potential and covert them to the correction field

mincy, cx = get_grad(corr_pot)
show_grad((mincy, cx), "\ncorrection pot-l")

# the correction is the exchanged components

show_grad((-cx[2:(end - 1), :], mincy[:, 2:(end - 1)]), "\ncorrector")

# Now we can subtract them
phix = g1 + cx[2:(end - 1), :]
phiy = g2 - mincy[:, 2:(end - 1)]
show_grad((phix, phiy), "\nCorrected phase")

# The corrected field should be integratable, what we can check by computing the cross derivatives
dxx, dxy = get_grad(phix)
dyx, dyy = get_grad(phiy)
fig, ax, hm = hmap(dxy - dyx);
ax.title = "mixed derivative after correction"
Colorbar(fig[1, 2], hm)
fig |> display

# ### Integration of the consisten gradient field in the aperture
#
# We start by integrating the corrected gradients along the aperture contour
cont, apedge = find_cw_border(ap)
ci2Point(ci) = Point(ci.I)
edgepos = ci2Point.(cont)
edgedirs = ci2Point.(circshift(cont, -1) .- cont)

function integrate_along_path(path, gr, x0=0)
    ret = zeros(length(path))
    ret[1] = x0
    dirs = circshift(path, -1) .- path
    for i in 1:(length(path) - 1)
        grad_ind = minimum(path[i:(i + 1)])
        grad = getindex.(gr, Ref(grad_ind))
        dir = [dirs[i].I...]
        important = dir .!= 0
        delta = dir[important]' * grad[important] # to avoid 0 * NaN
        ret[i + 1] = ret[i] + delta
    end
    return ret
end

function integrate_along_path_cyclic(path, gr)
    dirs = path .- circshift(path, 1)
    deltas = zeros(length(path))
    for i in 1:(length(path))
        j = mod1(i + 1, length(path))
        grad_ind = minimum(path[[i, j]])
        grad = getindex.(gr, Ref(grad_ind))
        dir = [dirs[i].I...]
        important = dir .!= 0
        deltas[i] = dir[important]' * grad[important] # to avoid 0 * NaN
    end
    ret = integrate_periodic_grad(deltas)
    return ret
end


edgeval = integrate_along_path(cont, [phix, phiy], ph[cont[1]])

lims = extrema(edgeval)

fig, ax, hm = hmap(phm; colorrange=lims);
arrows!(edgepos, 0.5edgedirs; color=:white)
scatterlines!(edgepos; color=:white)
ax.title = "GT phase and aperture contour"
Colorbar(fig[1, 2], hm)
fig |> display

scatter!(edgepos; color=edgeval, markersize=7.5, colorrange=lims)
ax.title = "GT phase and unwrapped phase on the aperture contour"
fig |> display

# Now let's replace in the gradients all element corresponding to the "left" edge (that is `[0,-1]` orientation) with the values of the restored aperture
scatter!(ci2Point.(apedge[:left]) .- Ref([0, 0.3]); marker='↓', color=:red)
fig |> display

# These points we want to update
ind = apedge[:left] .+ Ref(CartesianIndex(0, -1))

fig, ax, hm = hmap(phiy)
scatter!(ci2Point.(ind); marker='O', color=:red)
fig |> display

# Create dictionary to make replacement easier
#
edgedict = Dict(zip(cont, edgeval))

phiy[ind] .= getindex.(Ref(edgedict), apedge[:left])

# Repeat for the right edge, but with sign -1
# # These points we want to update
ind = apedge[:right] .+ Ref(CartesianIndex(0, 0))

fig, ax, hm = hmap(phiy)
scatter!(ci2Point.(ind); marker='O', color=:red)
fig |> display


phiy[ind] .= -getindex.(Ref(edgedict), apedge[:right])

# And for the up and down edges
ind = apedge[:up] .+ Ref(CartesianIndex(-1, 0))

fig, ax, hm = hmap(phix)
scatter!(ci2Point.(ind); marker='O', color=:red)
fig |> display


phix[ind] .= getindex.(Ref(edgedict), apedge[:up])

ind = apedge[:down]

fig, ax, hm = hmap(phix)
scatter!(ci2Point.(ind); marker='O', color=:red)
fig |> display


phix[ind] .= -getindex.(Ref(edgedict), apedge[:down])

# Replace the rest with zeroes, and we can integrate
phix[isnan.(phix)] .= 0
phiy[isnan.(phiy)] .= 0
show_grad((phix, phiy), "\nCorrected for integration")

# We'll inegrate using Finite difference cyclic, for this we need to append (prepend?) a raw and a column of zeroes
# phixcirc = vcat(zeros(size(ph, 1))', phix)
# phiycirc = hcat(zeros(size(ph, 1)), phiy)
phixcirc = vcat(phix, zeros(size(ph, 2))')
phiycirc = hcat(phiy, zeros(size(ph, 1)))

show_grad((phixcirc, phiycirc), "\nCorrected for integration")

# ## Integration by summation
#
# Now, we can intagrate simply by summation.
# The corrected gradient field should have zero residues.


gradres = PhaseUtils.getresmap(phixcirc, phiycirc)
fig, ax, hm = hmap(gradres);
ax.title = "Divergence of the corrected field"
Colorbar(fig[1, 2], hm)
fig |> display

ph_sum = circshift(reshape(cumsum(phixcirc[:]), size(ph)), 1)
fig, ax, hm = hmap(ph_sum)
ax.title = "Restored  by direct summation phase"
Colorbar(fig[1, 2], hm)
fig |> display

#  And we can calcualte the rotational component

rot_comp = phwrap(phm - ph_sum)
fig, ax, hm = hmap(rot_comp .* mask)
ax.title = "Rotational component"
Colorbar(fig[1, 2], hm)
fig |> display



# Restore with LS integration.

phi = circshift(
    integrate_2dgrad(phixcirc, phiycirc, PhaseUtils.FiniteDifferencesCyclic()), (0, 0)
)
fig, ax, hm = hmap(phi)
ax.title = "Restored phase with LS method"
Colorbar(fig[1, 2], hm)
fig |> display

# Now the region outside the boundary is  flat
fig, ax, hm = hmap(phi .* mask)
ax.title = "Restored phase with LS method"
Colorbar(fig[1, 2], hm)
fig |> display


# This method is now realized as `_unwrap_LS_Poisson` and is called by default by `unwrap_LS`
phm_un = unwrap_LS(phm, ap; restore_piston=true)
fig, ax, hm = hmap(phm_un .* mask);
ax.title = "Restored phase with PhaseUtils method"
Colorbar(fig[1, 2], hm)
fig |> display
