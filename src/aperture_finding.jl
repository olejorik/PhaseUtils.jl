# Algorithms to detect aperture from one or several images



# 1. Aperture from several interferograms

find_aperture(imgs, alg::AbstractApertureDetectionAlg) =
    error("method of finding aperture with $(typeof(alg)) is not defined yet")

struct MeanIntensity{T} <: AbstractApertureDetectionAlg where {T<:AbstractThresholdAlg}
    threshold::T
end

MeanIntensity(th::Real) = MeanIntensity(HardThreshold(th, true)) # by default use relative hard threshold

function find_aperture(imgs, alg::MeanIntensity)

    # img = mean(x -> Gray.(x), [i for i in igrams])
    img = mean(imgs)
    m = zeros(Bool, size(img))
    th = alg.threshold(img)
    m[img .> th] .= true

    return m
end


struct VarIntensity{T} <: AbstractApertureDetectionAlg where {T<:AbstractThresholdAlg}
    threshold::T
end

VarIntensity(th::Real) = VarIntensity(HardThreshold(th, true))

function find_aperture(imgs, alg::VarIntensity)

    # img = mean(x -> Gray.(x), [i for i in igrams])
    img = var(imgs)
    m = zeros(Bool, size(img))
    th = alg.threshold(img)
    m[img .> th] .= true

    return m
end


struct MaxIntensity{T} <: AbstractApertureDetectionAlg where {T<:AbstractThresholdAlg}
    threshold::T
end

MaxIntensity(th::Real) = MaxIntensity(HardThreshold(th, true))

function find_aperture(imgs, alg::MaxIntensity)

    # img = mean(x -> Gray.(x), [i for i in igrams])
    img = copy(imgs[1])
    for ind in eachindex(img)
        img[ind] = maximum(getindex.(imgs, ind))
    end
    m = zeros(Bool, size(img))
    th = alg.threshold(img)
    m[img .> th] .= true

    return m
end
