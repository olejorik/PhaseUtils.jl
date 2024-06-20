# Algorithms to detect aperture from one or several images



# 1. Aperture from several interferograms

find_aperture(imgs, alg::AbstractApertureDetectionAlg) =
    error("method of finding aperture with $(typeof(alg)) is not defined yet")

struct MeanIntensity{T::AbstractThresholdAlg} <: AbstractApertureDetectionAlg
    threshold::T
end
