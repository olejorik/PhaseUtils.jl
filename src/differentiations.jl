struct Fourier end
struct FiniteDifferences end
struct FiniteDifferencesCyclic end

_calculate_gradient(a) = _calculate_gradient(a, FiniteDifferencesCyclic())

"""
    _calculate_gradient(a [, method] )

Document this function
"""
function _calculate_gradient(a, ::Fourier)
    return real([ifft(fft(a, i) .* _grad_kernel(a, i), i) for i in 1:ndims(a)])
end

"""
    _gradfreq(n)

Calculate frequencies suitable for gradient calculations 
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

_grad_kernel(a, i) = _grad_kernel(a, i, Fourier())

function _grad_kernel(a, i, ::FiniteDifferencesCyclic)
    _dim = ones(Int, ndims(a))
    _dim[i] = size(a)[i]
    return 1 .- reshape(exp.(-2π .* im .* fftfreq(size(a)[i])), Tuple(_dim))
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
    _calculate_Laplacian(a, ::FiniteDifferencesCyclic)

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

_calculate_Laplacian(a) = _calculate_Laplacian(a, FiniteDifferencesCyclic())
