# Phase unwrapping algorithms

export itoh, unwrap_LS

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

Unwrap 2D `phase` defined inside `aperture`` using the Least-Squares decomposition
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
        piston = sum(sol[findall(ap .== 0)]) / length(findall(ap .== 0))
        sol .-= piston # value outside the aperture should be 0
    end

    return sol
end
