abstract type GradientCalculationAlg end

struct Fourier <: GradientCalculationAlg end
struct FiniteDifferences <: GradientCalculationAlg end
struct FiniteDifferencesCyclic <: GradientCalculationAlg end

default_grad_method(a) = FiniteDifferencesCyclic()

"""
    _calculate_gradient(a [, method] )

Return components of the gradient of the array`a` calculated using
`method`(`FiniteDifferencesCyclic` by default).
"""
_calculate_gradient(a) = _calculate_gradient(a, default_grad_method(a))

function _calculate_gradient(a, method::GradientCalculationAlg)
    return error("Implement $(typeof(method)) for gradient calcualtion")
end

function _calculate_gradient(a, ::Fourier)
    return real([ifft(fft(a, i) .* _grad_kernel(a, i), i) for i in 1:ndims(a)])
end

function _calculate_gradient(a, ::FiniteDifferences)
    return [a[2:end, :] .- a[1:(end - 1), :], a[:, 2:end] .- a[:, 1:(end - 1)]]
end

function _calculate_gradient(a, ::FiniteDifferencesCyclic)
    # ax = copy(a)
    # ay = copy(a)

    # ax[2:end, :] .= a[2:end, :] .- a[1:(end - 1), :]
    # ax[1, :] = a[1, :] .- a[end, :]
    # ay[:, 2:end] .= a[:, 2:end] .- a[:, 1:(end - 1)]
    # ay[:, 1] .= a[:, 1] .- a[:, end]

    ax = a .- circshift(a, Tuple(-step.left))
    ay = a .- circshift(a, Tuple(-step.up))
    return [ax, ay]
end

"""
    _gradfreq(n)

Calculate frequencies suitable for gradient calculations. These are the
    `fftfreq`, with the largest frequency (the Nyquist frequency) set to zero
    so the Fourier transform of the gradient calculated from a real function
    is Hermitian.
"""
function _gradfreq(n)
    f = collect(fftfreq(n))
    if iseven(n)
        f[n ÷ 2 + 1] = 0
    end
    return f
end  # function _gradfreq

function _grad_kernel(a, i, ::Fourier)
    _dim = ones(Int, ndims(a))
    _dim[i] = size(a)[i]
    return reshape(-im .* _gradfreq(size(a)[i]), Tuple(_dim))
end

function _grad_kernel(a, i, ::FiniteDifferencesCyclic)
    _dim = ones(Int, ndims(a))
    _dim[i] = size(a)[i]
    return 1 .- reshape(exp.(-2π .* im .* fftfreq(size(a)[i])), Tuple(_dim))
end

_grad_kernel(a, i) = _grad_kernel(a, i, default_grad_method(a))

"""
    _calculate_Laplacian(a [, method::GradientCalculationAlg])

TODO: add other methods and description
"""
function _calculate_Laplacian(a, ::FiniteDifferencesCyclic)
    ax2 = copy(a)
    ay2 = copy(a)

    ax2[2:(end - 1), :] .= 2a[2:(end - 1), :] .- a[1:(end - 2), :] .- a[3:end, :]
    ax2[1, :] = 2a[1, :] .- a[end, :] .- a[2, :]
    ax2[end, :] = 2a[end, :] .- a[end - 1, :] .- a[1, :]
    ay2[:, 2:(end - 1)] .= 2a[:, 2:(end - 1)] .- a[:, 1:(end - 2)] .- a[:, 3:end]
    ay2[:, 1] .= 2a[:, 1] .- a[:, end] .- a[:, 2]
    ay2[:, end] .= 2a[:, end] .- a[:, end - 1] .- a[:, 1]
    return ax2 .+ ay2
end

function _calculate_Laplacian(a, method::GradientCalculationAlg)
    return error("Implement $(typeof(method)) for Laplacian calculation")
end

_calculate_Laplacian(a) = _calculate_Laplacian(a, default_grad_method(a))
