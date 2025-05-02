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

    # ret = zeros(promote_type(eltype(vals), typeof(defaultel)), retsize)
    # if defaultel != 0
    #     ret .= defaultel
    # end
    # offset = CartesianIndex(to_indices(ret, 1 .- startoffset .+ shift))
    # ret[ind .+ offset] .= vals

    # return ret

    return _toArray(ind, vals, retsize, startoffset; shift=shift, defaultel=defaultel)

end

toArray(ind::Array{Tuple}, args...; kwargs...) =
    toArray(CartesianIndex.(ind), args...; kwargs...)

toArray(ind::Array{T}, vals, retsize::Tuple; kwargs...) where {T<:CartesianIndex} =
    _toArray(ind, vals, retsize, Tuple(oneunit(T)); kwargs...)


function _toArray(
    ind::Array{T}, vals, retsize::Tuple, startoffset::Tuple; shift=0, defaultel=0
) where {T<:CartesianIndex}
    ret = zeros(promote_type(eltype(vals), typeof(defaultel)), retsize)
    if defaultel != 0
        ret .= defaultel
    end
    offset = CartesianIndex(to_indices(ret, 1 .- startoffset .+ shift))
    ret[ind .+ offset] .= vals

    return ret
end
