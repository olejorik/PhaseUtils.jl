# Utils to vectorize and devectroize the data


"""
    toArray(ind<:Array{CartesianIndex}, vals; crop = true, shift=0, defaultel = 0) --> A

Construct an array `A` for which A[ind[i]] = val[i]. If `crop` is true, return only the minimal bounding box of `A`. shifted to `shift` pixels.
"""
function toArray(
    ind::Array{T}, vals; crop=true, shift=0, defaultel=0
) where {T<:CartesianIndex}

    ex = extrema(ind)
    if crop
        startoffset = Tuple(ex[1])
    else
        startoffset = Tuple(oneunit(T))
    end

    retsize = Tuple(ex[2]) .- startoffset .+ shift .+ 1

    ret = zeros(promote_type(eltype(vals), typeof(defaultel)), retsize)
    if defaultel != 0
        ret .= defaultel
    end
    offset = CartesianIndex(to_indices(ret, 1 .- startoffset .+ shift))
    ret[ind .+ offset] .= vals

    return ret


end

toArray(ind::Array{Tuple}, vals; kwargs...) = toArray(CartesianIndex.(ind), vals; kwargs...)
