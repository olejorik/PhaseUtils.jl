
function integrate_2dgrad(gx, gy)
    GX = fft(gx)
    GY = fft(gy)

    kxx = _grad_kernel(gx, 2, FiniteDifferencesCyclic())
    kyy = _grad_kernel(gy, 1, FiniteDifferencesCyclic())

    solfft = (GX .* conj(kxx) .+ GY .* conj(kyy)) ./ (abs2.(kxx) .+ abs2.(kyy))
    solfft[1, 1] = 0

    return real(ifft(solfft))
end

export integrate_2dgrad

