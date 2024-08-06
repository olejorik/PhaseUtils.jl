module SimpleProjections
import ..dotproduct
using LinearAlgebra: dot

export project, project!, reflect, reflect!, angle

function project!(a::AbstractArray, b::AbstractArray)
    @assert size(a) == size(b) "Arrays should have the same size"

    a[:] = dot(b, a) / dot(b, b) .* b[:]
    return a
end

project(a::AbstractArray, b::AbstractArray) = project!(copy(a), b)

function reflect!(a, b)
    p = project(a, b)
    a[:] .= 2 .* p[:] .- a[:]
    return a
end

reflect(a, b) = reflect!(copy(a), b)

angle(a, b) = acos(dot(a, b) / sqrt(dot(a, a) * dot(b, b)))




end
