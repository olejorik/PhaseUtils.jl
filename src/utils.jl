
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
diffirst(v, i) = [t .- v[i] for t in v[vcat(1:(i - 1), (i + 1):end)]]

diffirst!(target, v) = begin
    for i in eachindex(target)
        target[i] .= v[i + 1] .- v[1]
    end
    return target
end
diffirst!(target, v, i) = begin
    idxs = vcat(1:(i - 1), (i + 1):length(v))
    for (j, idx) in enumerate(idxs)
        target[j] .= v[idx] .- v[i]
    end
    return target
end

dotproduct(a, b) = sum([xa .* xb for (xa, xb) in zip(a, b)])

maskedrmse(a, binmask) = sqrt(sum(abs2, a[Bool.(binmask)]) / sum(binmask))
maskedrmse(a, b, binmask) = maskedrmse(a .- b, binmask)

maskedphasermse(a, binmask) = maskedrmse(phwrap(a), binmask)
maskedphasermse(a, b, binmask) = maskedphasermse(a .- b, binmask)

maskedPV(a, binmask) = (x -> last(x) - first(x))(extrema(a[Bool.(binmask)]))

"""
    ap2mask(ap)

Converts  `ap` array to mask, so that every zero is mapped to `NaN` and non-zero elements to 1.

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

"""
    binarize(ap)

Converts  `ap` array to a binary array `b`, so that every zero or `NaN` is mapped to 0 and non-zero elements to 1.
"""
function binarize(ap)
    b = ones(Bool, size(ap))
    b[isnan.(ap)] .= 0
    b[ap .== 0] .= 0
    return b
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
    circlemask!(a::Matrix{<:Number}, r, policy::ArrayAxes)

Zero elements of `a` outside the circle of radius `r` using coordinates from
`policy(size(a))`. The circle is centered at `0` if that value exists in each
axis; otherwise at the middle coordinate value.
"""
function circlemask!(a::Matrix{<:Number}, r, policy::ArrayAxes)
    axes = policy(size(a))
    cx = any(==(0), axes[1]) ? 0 : axes[1][cld(length(axes[1]), 2)]
    cy = any(==(0), axes[2]) ? 0 : axes[2][cld(length(axes[2]), 2)]
    return circlemask!(a, axes, cx, cy, r)
end

"""
    circlemask!(a::Matrix{<:Number}, cx, cy, r, policy::ArrayAxes)

Zero elements of `a` outside the circle centered at `(cx, cy)` with radius `r`
using coordinates from `policy(size(a))`.
"""
function circlemask!(a::Matrix{<:Number}, cx, cy, r, policy::ArrayAxes)
    axes = policy(size(a))
    return circlemask!(a, axes, cx, cy, r)
end

"""
    circlemask!(a::Matrix{<:Number}, axes::Vector{<:AbstractVector}, cx, cy, r)

Zero elements of `a` outside the circle centered at `(cx, cy)` with radius `r`
evaluated on explicit `axes` (one vector per dimension). The length of each
axis must match the corresponding dimension of `a`.
"""
function circlemask!(
    a::Matrix{<:Number}, axes::Vector{T}, cx, cy, r
) where {T<:AbstractVector}
    @assert length(axes) == 2 "Provide two axis vectors for a 2D array"
    @assert length(axes[1]) == size(a, 1) && length(axes[2]) == size(a, 2) "Axes must match array size"
    z = zero(eltype(a))
    for i1 in eachindex(axes[1]), i2 in eachindex(axes[2])
        ((axes[1][i1] - cx)^2 + (axes[2][i2] - cy)^2 > r^2) && (a[i1, i2] = z)
    end
    return nothing
end

"""
    circlemask(dims::NTuple{2,Int}, r, policy::ArrayAxes)

Create a boolean mask of size `dims` using coordinates provided by `policy`.
The circle is centered at the default origin of those coordinates (zero if
present; otherwise the middle coordinate value) with radius `r` expressed in
the same coordinate units.
"""
function circlemask(dims::NTuple{2,Int}, r, policy::ArrayAxes)
    axes = policy(dims)
    cx = any(==(0), axes[1]) ? 0 : axes[1][cld(length(axes[1]), 2)]
    cy = any(==(0), axes[2]) ? 0 : axes[2][cld(length(axes[2]), 2)]
    return [((x - cx)^2 + (y - cy)^2) <= r^2 for x in axes[1], y in axes[2]]
end

"""
    circlemask(axes::Vector{<:AbstractVector}, cx, cy, r)

Create a boolean mask evaluated on explicit coordinate `axes` (one vector per
dimension), marking points inside the circle of center `(cx, cy)` and radius
`r`.
"""
function circlemask(axes::Vector{T}, cx, cy, r) where {T<:AbstractVector}
    @assert length(axes) == 2 "Provide two axis vectors for a 2D mask"
    return [((x - cx)^2 + (y - cy)^2) <= r^2 for x in axes[1], y in axes[2]]
end

"""
    circlemask(dims::NTuple{2,Int}, cx, cy, r, policy::ArrayAxes)

Create a boolean mask of size `dims` using coordinates from `policy`, with a
circle centered at `(cx, cy)` in those coordinate units and radius `r`.
"""
function circlemask(dims::NTuple{2,Int}, cx, cy, r, policy::ArrayAxes)
    axes = policy(dims)
    return [((x - cx)^2 + (y - cy)^2) <= r^2 for x in axes[1], y in axes[2]]
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

"""
    linearphase(dims::NTuple{2,Int}, t::Tilt, policy::ArrayAxes=FourierAxes())

Create a matrix by evaluating the tilt `t` on coordinates given by `policy(dims)`.
This is equivalent to `materialize(t, policy(dims))`.
"""
function linearphase(dims::NTuple{2,Int}, t::Tilt, policy::ArrayAxes=FourierAxes())
    return materialize(t, policy(dims))
end

"""
    linearphase(T::Type{<:Number}, dims::NTuple{2,Int}, t::Tilt, policy::ArrayAxes=FourierAxes())

Typed variant of `linearphase(dims, t, policy)` converting the result to `T`.
"""
function linearphase(
    T::Type{<:Number}, dims::NTuple{2,Int}, t::Tilt, policy::ArrayAxes=FourierAxes()
)
    return convert.(T, materialize(t, policy(dims)))
end

"""
    linearphase(axes::Vector{<:AbstractVector}, t::Tilt)

Evaluate the tilt `t` on explicit coordinate `axes` (one vector per dimension).
"""
function linearphase(axes::Vector{T}, t::Tilt) where {T<:AbstractVector}
    return materialize(t, axes)
end
