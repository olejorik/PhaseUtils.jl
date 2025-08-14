module PhaseUtils
using FFTW
using FileIO
using Statistics

export phwrap, maskedrmse, maskedphasermse, ap2mask, mask2ap, binarize
export bboxview
export hardthreshold, hardthreshold!, softthreshold
export circlemask, circlemask!, linearphase
export toArray
export readfiles
export find_aperture
export findpiston, findtiptilt
export Tilt,
    TiltCentered, FreeTilt, sigma, tau, setsigma!, settau!, setall!, materialize, apply
export ArrayAxes, FourierAxes, DataAxes, DataAxesCentered, ProvidedAxes



include("axes_and_tilts.jl")
include("utils.jl")
include("fileutils.jl")

include("algorithms.jl")
include("aperture_border.jl")
include("aperture_finding.jl")
include("differentiations.jl")
include("integrations.jl")
include("unwrapping.jl")
include("cropandpad.jl")
include("imageprocessing.jl")
include("vectorization.jl")
include("SimpleProjections.jl")
include("PTT.jl")

end
