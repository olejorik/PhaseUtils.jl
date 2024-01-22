
function igram(phase, a=0.5, b=0.5)
    return a .+ b .* cos.(phase)
end

"""
    phwrap(ϕ)

Wrap phase `ϕ`: ψ = `phwrap`(ϕ) ⇔ ψ = ϕ + 2πk, k∈Z, -π <ψ ≤π.
If called on an array, works pointwise.
"""
phwrap(x::Number) = isnan(x) ? NaN : rem2pi(x, RoundNearest)

phwrap(a::AbstractArray) = phwrap.(a)
phwrap(::Missing) = missing

diffirst(v) = [t .- v[1] for t in v[2:end]]
dotproduct(a, b) = sum([xa .* xb for (xa, xb) in zip(a, b)])

maskedrmse(a, binmask) = sqrt(sum(abs2, a .* binmask) / sum(binmask))
maskedrmse(a, b, binmask) = maskedrmse(a .- b, binmask)

function maskedphasermse(a, b, binmask)
    return sqrt(sum(abs2, phwrap.((a .- b) .* binmask)) / sum(binmask))
end

"""
    ap2mask(ap)

Converts  `ap` array to mask, so that every zrro is mapped to `NaN` and non-zero elements to 1.

See also [`mask2ap`](@ref).
"""
function ap2mask(ap)
    # mask = Array{Union{Missing, Int}}(missing, size(ap)...)
    mask = Array{Float64}(undef, size(ap)...)
    mask[ap .!= 0] .= 1
    mask[ap .== 0] .= NaN
    return mask
end

"""
    mask2ap(mask)

Converts  `mask` array where points outside the aperture are defines as `NaN` to a Float array, so that `NaN` -> 0, `not NaN` -> 1.
"""
function mask2ap(mask)
    ap = zero(mask) .+ 1
    ap[isnan.(mask)] .= 0
    return ap
end

function bboxview(arr)
    idx = findall(!isnan, arr)
    region = minimum(idx):maximum(idx)
    return @view arr[region]
end

"""
    bboxview(arr, mask, pad = 0)

Make a box corresponding to the mask not-NaN elements surrounded by pad.
"""
function bboxview(arr, mask, pad=0)
    # idx = findall(!ismissing, mask)
    idx = findall(!isnan, mask)
    cpad = CartesianIndex(pad, pad)
    region = (minimum(idx) - cpad):(maximum(idx) + cpad)
    return @view arr[region]
end



# Small convenience utils
zerostomissing(a) = replace(a, 0 => missing)
missingtozeros(a) = replace(a, missing => 0)
equalizemissing(a) = replace(x -> ismissing(x) ? missing : 1, a)
missingtobinary(a) = missingtozeros(equalizemissing(a))

"""
    circlemask!(a::Matrix{<:Number}, [cx, cy,] r)

Set to zero values of matrix `a` outside the circle with center (`cx`,`cy`) and radius `r`.

If omitted, (`cx`,`cy`) is the center of the matrix.
"""
function circlemask!(a::Matrix{<:Number}, cx, cy, r)
    for i in eachindex(IndexCartesian(), a)
        ((i[1] - cx)^2 + (i[2] - cy)^2 > r^2) && (a[i] = 0)
    end
    return nothing
end

function circlemask!(a::Matrix{<:Number}, r)
    cx, cy = size(a) ./ 2 .+ 0.5
    return circlemask!(a, cx, cy, r)
end

"""
    circlemask(dims::NTuple{2, Int}, cx, cy, r)

Create a boolean matrix of size `dims` with `true` only inside the circle with center (`cx`,`cy`) and radius `r`.

If omitted, (`cx`,`cy`) is the center of the domain: `(dims+1)/2`.

## Examples
```julia
julia> circlemask((6,8), 2.5,3, 1.5)
6×8 Matrix{Bool}:
 0  0  1  0  0  0  0  0
 0  1  1  1  0  0  0  0
 0  1  1  1  0  0  0  0
 0  0  1  0  0  0  0  0
 0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0

julia> circlemask((10,10), 3)
10×10 Matrix{Bool}:
 0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0
 0  0  0  1  1  1  1  0  0  0
 0  0  1  1  1  1  1  1  0  0
 0  0  1  1  1  1  1  1  0  0
 0  0  1  1  1  1  1  1  0  0
 0  0  1  1  1  1  1  1  0  0
 0  0  0  1  1  1  1  0  0  0
 0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0
```
"""
function circlemask(dims::NTuple{2,Int}, args...)
    a = ones(Bool, dims)
    circlemask!(a, args...)
    return a
end

"""
    linearphase(T=Float64, dims::NTuple{2,Int},  cx, cy, kx, ky)

Create a matrix with linear phase with slopes (`kx`, `ky`) and having zero at point (`cx`, `cy`).

## Example
```julia
julia> linearphase((5,7), 3, 3, 1.5, 1.5)
5×7 Matrix{Float64}:
 -6.0  -4.5  -3.0  -1.5  0.0  1.5  3.0
 -4.5  -3.0  -1.5   0.0  1.5  3.0  4.5
 -3.0  -1.5   0.0   1.5  3.0  4.5  6.0
 -1.5   0.0   1.5   3.0  4.5  6.0  7.5
  0.0   1.5   3.0   4.5  6.0  7.5  9.0
```
"""
function linearphase(T::Type{<:Number}, dims::NTuple{2,Int}, cx, cy, kx, ky)
    a = Array{T,2}(undef, dims)
    for i in eachindex(IndexCartesian(), a)
        a[i] = (i[1] - cx) * kx + (i[2] - cy) * ky
    end
    return a
end

linearphase(dims::NTuple{2,Int}, cx, cy, kx, ky) =
    linearphase(Float64, dims, cx, cy, kx, ky)
