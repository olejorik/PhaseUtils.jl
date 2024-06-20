function hardthreshold!(a::AbstractArray, th)
    a[abs.(a) .<= th] .= 0
    return a
end

function hardthreshold!(a::AbstractArray, th, relative::Bool)
    if relative
        return hardthreshold!(a, th * maximum(a))
    else
        return hardthreshold!(a, th)
    end
end

hardthreshold(a, args...) = hardthreshold!(copy(a), args...)

_softthreshold(a, th) = sign.(a) .* max.(abs.(a) .- th, 0)

function softthreshold!(a::AbstractArray, th)
    a .= _softthreshold.(a, th)
    return a
end

function softthreshold!(a::AbstractArray, th, relative::Bool)
    if relative
        return softthreshold!(a, th * maximum(a))
    else
        return softthreshold!(a, th)
    end
end

function softthreshold(a::AbstractArray, th)
    return sign.(a) .* max.(abs.(a) .- th, 0)
end

function softthreshold(a::AbstractArray, th, relative::Bool)
    if relative
        return softthreshold(a, th * maximum(a))
    else
        return softthreshold(a, th)
    end
end


threshold!(a, alg::AbstractThresholdAlg) =
    error("method of thresholding with $(typeof(alg)) is not defined yet")

threshold(a, args...) = threshold!(copy(a), args...)

struct HardThreshold{T} <: AbstractThresholdAlg
    th::T
    relative::Bool
end

HardThreshold(th) = HardThreshold(th, false)

(alg::HardThreshold)(img) = alg.relative ? alg.th * maximum(img) : alg.th

threshold!(a, alg::HardThreshold) = hardthreshold!(a, alg.th, alg.relative)

struct SoftThreshold{T} <: AbstractThresholdAlg
    th::T
    relative::Bool
end

SoftThreshold(th) = SoftThreshold(th, false)

threshold!(a, alg::SoftThreshold) = softthreshold!(a, alg.th, alg.relative)
