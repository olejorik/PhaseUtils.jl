export crop

"""
    crop(a::AbstractArray, cropsize[, cropcenter])
    crop(cropsize[, cropcenter]) creates a cropping operator.

Crop the array to `cropsize` around `cropcenter`. If `cropcenter` is omitted, use central
element of the array.

`cropsize` and `cropcenter` can be either tuples/vectors of integers
or just single integer, which represent a tuple of identical numbers.

# Examples
```jldoctest
julia> a = reshape(1:15, 3, 5); crop(a, 1, (1, 3))
1×1 Matrix{Int64}:
 7

julia> crop(a, (2, 2), (1, 3))
1×2 Matrix{Int64}:
 4  7

julia> crop(3)(a)
3×3 Matrix{Int64}:
 4  7  10
 5  8  11
 6  9  12

```

"""
function crop(a::AbstractArray, cropsize, cropcenter)
    return a[intersect(CartesianIndices(a), croprange(Tuple(cropsize), cropcenter))]
end

crop(a::AbstractArray, cs::Int, cc) = crop(a, repeat([cs], ndims(a)), cc)

crop(a::AbstractArray, cropsize) = crop(a, cropsize, size(a) .÷ 2 .+ 1)
crop(cropsize, cropcenter) = (a -> crop(a, cropsize, cropcenter))
crop(cropsize) = (a -> crop(a, cropsize))

function croprange(cropsize::Tuple{Vararg{Int}}, cropcenter)
    mincorner = CartesianIndex((cropcenter .- cropsize .÷ 2)...)
    maxcorner = CartesianIndex((cropcenter .- cropsize .÷ 2 .+ cropsize .- 1)...)
    return mincorner:maxcorner
end
