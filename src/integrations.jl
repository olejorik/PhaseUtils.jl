
"""
    integrate_2dgrad(gx, gy[, gradmethod=default_grad_method(gx)])

TBW
"""
function integrate_2dgrad(gx, gy, ::LeastSquares, gradmethod=default_grad_method(gx))
    GX = fft(gx)
    GY = fft(gy)

    kxx = _grad_kernel(gx, 2, gradmethod)
    kyy = _grad_kernel(gy, 1, gradmethod)

    solfft = (GX .* conj(kxx) .+ GY .* conj(kyy)) ./ (abs2.(kxx) .+ abs2.(kyy))
    solfft[1, 1] = 0

    return real(ifft(solfft))
end

function integrate_2dgrad(gx, gy, method::InverseProblemAlg)
    return error("Implement $(typeof(method)) for gradient integration")
end

integrate_2dgrad(gx, gy, args...) = integrate_2dgrad(gx, gy, LeastSquares(), args...)

"""
    integrate_periodic_grad(g)

Integrate 1D gradient with periodic boundary conditions
"""
function integrate_periodic_grad(g)
    G = fft(g)

    k = _grad_kernel(g, 1, FiniteDifferencesCyclic())

    solfft = (G .* conj(k)) ./ abs2.(k)
    solfft[1] = 0

    return real(ifft(solfft))
end

export integrate_2dgrad, integrate_periodic_grad
