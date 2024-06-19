# Utils to vectorize and devectroize the data


"""
    toArray(ind<:Array{CartesianIndex}, vals; crop = true, shift=0) --> A

Construct an array `A` for which A[ind[i]] = val[i]. If `crop` is true, return only the minimal bounding box of `A`. shifted to `shift` pixels.
"""
function toArray(ind::Array{T}, vals; crop=true, shift=0) where {T<:CartesianIndex}

    ex = extrema(ind)
    if crop
        startoffset = Tuple(ex[1])
    else
        startoffset = Tuple(one(T))
    end

    retsize = Tuple(ex[2]) .- startoffset .+ shift .+ 1

    ret = zeros(eltype(vals), retsize)
    offset = CartesianIndex(to_indices(ret, 1 .- startoffset .+ shift))
    ret[ind .+ offset] .= vals

    return ret


end
