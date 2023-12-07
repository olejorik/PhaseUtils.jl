function hardthreshold!(a::AbstractArray, th)
    a[abs.(a) .<= th] .= 0
    return a
end

function hardthreshold!(a::AbstractArray, th, relative::Bool)
    if relative
        return hardthreshold!(a, th * maximum(a))
    else
        return hardthreshold!(a, th)
    end
end

hardthreshold(a, args...) = hardthreshold!(copy(a), args...)

function softthreshold(a::AbstractArray, th)
    return sign.(a) .* max.(abs.(a) .- th, 0)
end

function softthreshold(a::AbstractArray, th, relative::Bool)
    if relative
        return softthreshold(a, th * maximum(a))
    else
        return softthreshold(a, th)
    end
end
