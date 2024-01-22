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
function unwrap_LS(phase, ap; restore_piston=true)
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

function getresmap(phase, positions=true)
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
    return posy, posx, resmap
end

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

function inner_res(posx, posy, ap)
    x = Int(posx - 0.5)
    y = Int(posy - 0.5)
    return all(Bool[ap[x, y], ap[x + 1, y], ap[x, y + 1], ap[x + 1, y + 1]])
end

function ideal_residue(x, y, resx, resy)
    return atan(y - resy, x - resx)
end
