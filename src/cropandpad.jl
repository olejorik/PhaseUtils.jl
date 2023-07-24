export crop

function crop(a::AbstractArray, cropsize, cropcenter)
    return a[intersect(CartesianIndices(a), croprange(Tuple(cropsize), cropcenter))]
end

crop(a, cropsize) = crop(a, cropsize, size(a) .÷ 2 .+ 1)

function croprange(cropsize::Tuple{Vararg{Int}}, cropcenter)
    mincorner = CartesianIndex((cropcenter .- cropsize .÷ 2)...)
    maxcorner = CartesianIndex((cropcenter .- cropsize .÷ 2 .+ cropsize .- 1)...)
    return mincorner:maxcorner
end
