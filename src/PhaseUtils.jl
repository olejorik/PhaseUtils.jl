module PhaseUtils
using FFTW

function igram(phase, a=0.5, b=0.5)
    return a .+ b .* cos.(phase)
end

"Wrap Phase"
phwrap(x::Number) = isnan(x) ? NaN : rem2pi(x, RoundNearest)

phwrap(a::AbstractArray) = phwrap.(a)
phwrap(::Missing) = missing

diffirst(v) = [t .- v[1] for t in v[2:end]]
dotproduct(a, b) = sum([xa .* xb for (xa, xb) in zip(a, b)])

maskedrmse(a, b, binmask) = sqrt(sum(abs2, (a .- b) .* binmask) / sum(binmask))
function maskedphasermse(a, b, binmask)
    return sqrt(sum(abs2, phwrap.((a .- b) .* binmask)) / sum(binmask))
end

function ap2mask(ap)
    mask = zero(ap) .+ 1
    mask[ap .== 0] .= missing
    return mask
end

function mask2ap(mask)
    ap = zero(mask) .+ 1
    ap[ismissing.(mask)] .= 0
    return ap
end

function bboxview(arr)
    idx = findall(!ismissing, arr)
    region = minimum(idx):maximum(idx)
    return @view arr[region]
end

"""
    bboxview(arr, mask, pad = 0)

Make a box corresponding to the mask not-NaN elements surrounded by pad.
"""
function bboxview(arr, mask, pad=0)
    idx = findall(!ismissing, mask)
    cpad = CartesianIndex(pad, pad)
    region = (minimum(idx) - cpad):(maximum(idx) + cpad)
    return @view arr[region]
end

export phwrap, maskedrmse, maskedphasermse, ap2mask, mask2ap
export bboxview

# Small convenience utils
zerostomissing(a) = replace(a, 0 => missing)
missingtozeros(a) = replace(a, missing => 0)
equalizemissing(a) = replace(x -> ismissing(x) ? missing : 1, a)
missingtobinary(a) = missingtozeros(equalizemissing(a))

include("aperture_border.jl")
include("differentiations.jl")
include("integrations.jl")
include("unwrapping.jl")

end
