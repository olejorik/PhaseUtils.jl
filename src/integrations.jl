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