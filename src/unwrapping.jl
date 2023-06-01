"""
    itoh(phi)

Itoh's algorithm for 1D phase unwrapping.
"""
function itoh(phi)
    dphi = phwrap(circshift(phi, -1) .- phi)
    phi_restored = similar(dphi)
    integral = phi[1]
    for i in eachindex(dphi)
        phi_restored[i] = integral
        integral += dphi[i]
    end
    return phi_restored
end  # function itohphi

export itoh

function unwrap_LS(phase, ap; restore_piston = true)
    # Calculate wrapped gradients
    kxx = _grad_kernel(phase, 1, FiniteDifferencesCyclic())
    kyy = _grad_kernel(phase, 2, FiniteDifferencesCyclic())

    phase_0 = phase .* ap
    wgx, wgy = _calculate_gradient(phase_0) |> phwrap


    # Calculate the boundary
    cw_cont_ap, edge =find_cw_border(ap)

    # unwrap the phase on the boundary
    ph_b = phase_0[cw_cont_ap]
    ph_r = itoh(ph_b)

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
        piston = sum(sol[findall(ap .== 0)])/length(findall(ap .== 0))
        sol .-= piston # value outside the aperture should be 0
    end
    
    return sol

end

export unwrap_LS
