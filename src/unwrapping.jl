# Phase unwrapping algorithms

export itoh, unwrap_LS, getresmap, getresmapsparce

"""
    itoh(phi)

Itoh's algorithm for 1D phase unwrapping.

# Example
```jldoctest
julia> ph = 1:10
1:10

julia> ph_w = phwrap(ph)
10-element Vector{Float64}:
  1.0
  2.0
  3.0
 -2.2831853071795867
 -1.2831853071795865
 -0.28318530717958645
  0.7168146928204135
  1.7168146928204135
  2.7168146928204133
 -2.566370614359173

julia> itoh(ph_w)
10-element Vector{Float64}:
  1.0
  2.0
  3.0
  4.0
  5.0
  6.0
  7.0
  8.0
  9.0
 10.0


```
"""
function itoh(phi)
    return _itoh(phi)[1]
end

"""
    _itoh(phi) -> itoh(phi), nres

Return unwrapped phase and sum of the residues inside the path if `phi` are the values
on its boundary.
"""
function _itoh(phi)
    dphi = phwrap(circshift(phi, -1) .- phi)
    phi_restored = similar(dphi)
    integral = phi[1]
    for i in eachindex(dphi)
        phi_restored[i] = integral
        integral += dphi[i]
    end
    return phi_restored, (phi[1] + integral) / 2π #number of the residues in the closed path
end

"""
    unwrap_LS(phase, aperture; restore_piston=true)

Unwrap 2D `phase` defined inside `aperture` using the Least-Squares decomposition
of the wrapped gradient of the wprapped phase in the rotor-free and solenodial field
and integration of the rotor-free part.

"""
unwrap_LS(phase, ap; restore_piston=false) = _unwrap_LS_Poisson(phase, ap; restore_piston)

function _unwrap_LS_quick(phase, ap; restore_piston=true)
    # Calculate wrapped gradients
    kxx = _grad_kernel(phase, 1, FiniteDifferencesCyclic())
    kyy = _grad_kernel(phase, 2, FiniteDifferencesCyclic())

    phase_0 = phase .* ap
    wgx, wgy = phwrap(_calculate_gradient(phase_0))

    # Calculate the boundary
    cw_cont_ap, edge = find_cw_border(ap)

    # unwrap the phase on the boundary
    ph_b = phase_0[cw_cont_ap]
    # ph_r = itoh(ph_b) # this is valid only for consistent case
    # in a general case, use the least-squares approach
    # ph_r = unwrap_contour_LS(ph_b)
    #
    # Unwrap phase on the boundary with correction for the residues
    #
    ph_r = copy(ph_b) #
    posx, posy, resmap = getresmapsparce(phase_0) # get all residues
    # @info "$(length(resmap)) residues found)"
    inrescount = 0
    for (px, py, r) in zip(posx, posy, resmap)
        if inner_res(px, py, ap)
            correction(ind) = ideal_residue(Tuple(ind)..., px, py)
            ph_r .-= r * correction.(cw_cont_ap)
            inrescount += 1
        end
    end
    # ph_r .= itoh(ph_r)
    ph_r = unwrap_contour_LS(ph_r)
    @info "$inrescount ineternal residues corrected for contour unwrapping"


    # Substite in the gradients the values calculated on the edges with the real ones

    edgedict = Dict(zip(cw_cont_ap, ph_r))

    for i in eachindex(edge.left)
        wgx[edge.left[i]] = edgedict[edge.left[i]]
    end

    for i in eachindex(edge.right)
        wgx[edge.right[i] + step.right] = -edgedict[edge.right[i]]
        # wgx[edge.right[i] ] = edgedict[edge.right[i]]
    end

    for i in eachindex(edge.up)
        wgy[edge.up[i]] = edgedict[edge.up[i]]
    end

    for i in eachindex(edge.down)
        wgy[edge.down[i] + step.down] = -edgedict[edge.down[i]]
        # wgy[edge.down[i] ] = edgedict[edge.down[i]]
    end

    # integrate the true gradients
    sol = integrate_2dgrad(wgx, wgy)

    if restore_piston
        piston = sum(sol[findall(ap .== 0)]) / length(findall(ap .== 0))
        sol .-= piston # value outside the aperture should be 0
    end

    return sol
end

function _unwrap_LS_Poisson(phase, ap; restore_piston=false)

    diffmeth = PhaseUtils.FiniteDifferences()
    get_grad(arr) = PhaseUtils._calculate_gradient(arr, diffmeth)

    # Calculate wrapped gradients

    mask = ap2mask(ap)
    phm = phase .* mask
    g = (g1, g2) = phwrap(get_grad(phm))
    posx, posy, resmap = getresmapsparce(phm)
    dualap_inside = dual_region_box(ap)
    rho = zeros(size(dualap_inside))
    for ind in eachindex(resmap)
        rho[dualcoordplus(posx[ind], posy[ind])] = π / 2 * resmap[ind]
    end

    corr_pot = membrane_sor(rho, dualap_inside)
    mincy, cx = get_grad(corr_pot)

    phix = g1 + cx[2:(end - 1), :]
    phiy = g2 - mincy[:, 2:(end - 1)]
    dxx, dxy = get_grad(phix)
    dyx, dyy = get_grad(phiy)


    # Calculate the boundary
    cont, apedge = find_cw_border(ap)
    # ci2Point(ci) = Point(ci.I)
    # edgepos = ci2Point.(cont)
    # edgedirs = ci2Point.(circshift(cont, -1) .- cont)

    # unwrap the phase on the boundary
    edgeval = integrate_along_path(cont, [phix, phiy], phm[cont[1]])


    # Unwrap phase on the boundary with correction for the residues
    #

    # Substite in the gradients the values calculated on the edges with the real ones
    # These points we want to update
    ind = apedge[:left] .+ Ref(CartesianIndex(0, -1))

    # Create dictionary to make replacement easier
    #
    edgedict = Dict(zip(cont, edgeval))

    phiy[ind] .= getindex.(Ref(edgedict), apedge[:left])

    # Repeat for the right edge, but with sign -1
    # # These points we want to update
    ind = apedge[:right] .+ Ref(CartesianIndex(0, 0))

    phiy[ind] .= -getindex.(Ref(edgedict), apedge[:right])

    # And for the up and down edges
    ind = apedge[:up] .+ Ref(CartesianIndex(-1, 0))

    phix[ind] .= getindex.(Ref(edgedict), apedge[:up])

    ind = apedge[:down]


    phix[ind] .= -getindex.(Ref(edgedict), apedge[:down])

    # Replace the rest with zeroes, and we can integrate
    phix[isnan.(phix)] .= 0
    phiy[isnan.(phiy)] .= 0

    # We'll inegrate using Finite difference cyclic, for this we need to append (prepend?) a raw and a column of zeroes
    # phixcirc = vcat(zeros(size(ph, 1))', phix)
    # phiycirc = hcat(zeros(size(ph, 1)), phiy)
    phixcirc = vcat(phix, zeros(size(phm, 1))')
    phiycirc = hcat(phiy, zeros(size(phm, 1)))


    # integrate the true gradients
    sol = integrate_2dgrad(phixcirc, phiycirc, PhaseUtils.FiniteDifferencesCyclic())

    if restore_piston
        piston = sum(sol[findall(ap .== 0)]) / length(findall(ap .== 0))
        sol .-= piston # value outside the aperture should be 0
    end

    return sol
end

function unwrap_contour_LS(phase; restore_piston=false)
    wg = phwrap(phase .- circshift(phase, 1))

    sol = integrate_periodic_grad(wg)

    if restore_piston
        piston = sum(phase) / length(phase)
        sol .+= piston # the mean of the restored phase is equal to the mean of the wrapped phase
    end

    return sol
    # return itoh(phase)

end

function resblock(A)
    return sum(
        phwrap, [A[2, 1] - A[1, 1], A[2, 2] - A[2, 1], A[1, 2] - A[2, 2], A[1, 1] - A[1, 2]]
    )
end

"""
    getresmap(phase) -> resmap

Find residues in the phase, which correspond to the residues in the wrapped phase gradients expressed in 2pi multiples.
"""
function getresmap(phase)
    I, J = size(phase)
    resmap = similar(phase, I - 1, J - 1)
    posx = similar(phase, I - 1, J - 1)
    posy = similar(phase, I - 1, J - 1)
    for i in 1:(I - 1), j in 1:(J - 1)
        resmap[i, j] = resblock(phase[i:(i + 1), j:(j + 1)]) ÷ (2π)
        posx[i, j] = i + 1 / 2
        posy[i, j] = j + 1 / 2
    end
    resmap[resmap .== 0] .= NaN
    return resmap
end

"""
    getresmap(gx, gy) -> resmap

Find residues in the phase, which correspond to the residues in the wrapped phase gradients expressed in 2pi multiples.
"""
function getresmap(gx, gy)
    I, J = size(gx, 2), size(gy, 1)
    resmap = similar(gx, I - 1, J - 1)
    posx = similar(gx, I - 1, J - 1)
    posy = similar(gx, I - 1, J - 1)
    for i in 1:(I - 1), j in 1:(J - 1)
        resmap[i, j] = gx[i, j] + gy[i + 1, j] - gx[i, j + 1] - gy[i, j]
        posx[i, j] = i + 1 / 2
        posy[i, j] = j + 1 / 2
    end
    # resmap[resmap .== 0] .= NaN
    return resmap
end

"""
    getresmapsparce(phase) -> posx, posy, resmap

Find residues in the phase, which correspond to the residues in the wrapped phase gradients expressed in 2pi multiples.
"""
function getresmapsparce(phase)
    I, J = size(phase)
    resmap = Int[]
    posx = Float16[]
    posy = Float16[]
    for i in 1:(I - 1), j in 1:(J - 1)
        res = div(resblock(phase[i:(i + 1), j:(j + 1)]), (2π), RoundNearest)
        if !isnan(res) && res != 0
            push!(resmap, res)
            push!(posx, i + 1 / 2)
            push!(posy, j + 1 / 2)
        end
    end
    return posx, posy, resmap
end

"""
    getresmapsparce(gx,gy) -> posx, posy, resmap

Find residues in the (inconsistent) gradient field.
"""
function getresmapsparce(gx, gy)
    I, J = size(gx, 2), size(gy, 1)
    resmap = Float64[]
    posx = Float16[]
    posy = Float16[]
    for i in 1:(I - 1), j in 1:(J - 1)
        res = gx[i, j] + gy[i + 1, j] - gx[i, j + 1] - gy[i, j]
        if !isnan(res) && !(1 + res ≈ 1) # res is approximately zero if compared with 1
            push!(resmap, res)
            push!(posx, i + 1 / 2)
            push!(posy, j + 1 / 2)
        end
    end
    return posx, posy, resmap
end

function inner_res(posx, posy, ap)
    x = Int(posx - 0.5)
    y = Int(posy - 0.5)
    return all(Bool[ap[x, y], ap[x + 1, y], ap[x, y + 1], ap[x + 1, y + 1]])
end

function ideal_residue(x, y, resx, resy)
    return atan(y - resy, x - resx)
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
