
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

"""
    membrane_sor!(u, f::Array{Real}, mask_internal::Array{Bool})

Solves Poisson equation Δ `u` = `f`, with boundary conditions set by `u`(x) for x ∉ `mask_internal` using Successive Overrelaxation method (see Numerical Recipes (www.cambridge.org/9780521880688), ch.20.5.1).
"""
function membrane_sor!(u, f::Array, mask_internal::Array{Bool,2}; maxits=1000, tol=1e-6)

    size(f) == size(mask_internal) || error("Size of the source and mask do not match")

    anormf = sum(abs, f) # Compute initial norm of residual and terminate iterations when norm has been reduced by a factor tol
    omega = 1.0
    jmax, imax = size(f)
    rjac = (cos(π / jmax) + cos(π / imax)) / 2
    anormlast = anormf

    for n in 1:maxits
        anorm = 0
        for ipass in (0, 1) #odd-even numbering
            lineshift = ipass
            for j in 2:(jmax - 1)
                for l in (lineshift + 2):2:(imax - 1)
                    if mask_internal[j, l]
                        # @info "j=$j, l=$l"
                        resid =
                            (u[j + 1, l] + u[j - 1, l] + u[j, l + 1] + u[j, l - 1]) / 4 -
                            u[j, l] - f[j, l]
                        anorm += abs(resid)
                        u[j, l] = u[j, l] + omega * resid
                    end
                end
                lineshift = 1 - lineshift
            end
            omega = (
                if n == 0 && ipass == 0
                    1.0 / (1.0 - 0.5rjac^2)
                else
                    1.0 / (1.0 - 0.25rjac^2 * omega)
                end
            )
        end

        anormlast = anorm

        if (anormlast < tol * anormf)
            @info "converged in $n iterations"
            @info "anormf = $anormf"
            @info "anormlast = $anormlast"
            @info "omega = $omega"
            return nothing
        end

    end
    @info "failed to converge in $maxits iterations"
    @info "anormf = $anormf"
    @info "anormlast = $anormlast"
    @info "omega = $omega"
    return nothing


end

"""
    membrane_sor(f::Array{Real}, mask_internal::Array{Bool})

Solves Poisson equation Δ `u` = `f`, with boundary conditions `u`(x) = 0 <=> x ∉ `mask_internal` using Successive Overrelaxation method (see Numerical Recipes (www.cambridge.org/9780521880688), ch.20.5.1).
"""
function membrane_sor(f::Array, mask_internal::Array{Bool,2}; kwargs...)
    u = zeros(size(f))
    membrane_sor!(u, f, mask_internal; kwargs...)
    return u
end

export integrate_2dgrad, integrate_periodic_grad, membrane_sor!, membrane_sor
