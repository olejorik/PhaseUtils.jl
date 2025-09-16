"""
    Windowing functionality for Phase project

This module provides windowing functions for image processing and signal analysis.
"""

abstract type Window{T<:Union{Real,Tuple}} end

(w::Window)(dims::Tuple{Vararg{Int}}) =
    error("Windowing for type $(typeof(w)) is not defined")

"""
    GaussianWindow(w)

Creates a centered Gaussian window `w` pixels wide.

`w` can be a single number or a tuple, with specified width for each dimension.

# Example
```julia
# Create a symmetric Gaussian window
gw = GaussianWindow(15)
window_2d = gw((64, 64))

# Create an asymmetric Gaussian window
gw_asym = GaussianWindow((20, 10))
window_asym = gw_asym((64, 32))
```
"""
struct GaussianWindow{T} <: Window{T}
    width::T
end

_grange(len, wid) = range(-(len - 1) / wid, (len - 1) / wid; length=len)
_gwidth(w::GaussianWindow{<:Real}, dims) = fill(w.width, length(dims))
_gwidth(w::GaussianWindow{<:Tuple}, dims) =
    if length(w.width) == length(dims)
        w.width
    else
        error("Incompatible dimensions of the window and the array")
    end

function (w::GaussianWindow)(dims::Tuple{Vararg{Int}})
    a = zeros(dims) .+ 1
    widths = _gwidth(w, dims)
    for (id, d) in enumerate(dims)
        wd = _grange(d, widths[id])
        for (is, s) in enumerate(eachslice(a; dims=id))
            s .*= exp.(-wd[is]^2 / 2)
        end
    end
    return a
end
