"""
    Windowing functionality for Phase project

This module provides windowing functions for image processing and signal analysis.
"""

abstract type Window{T<:Union{Real,Tuple}} end

(w::Window)(dims::Tuple{Vararg{Int}}) =
    error("Windowing for type $(typeof(w)) is not defined")

"""
    GaussianWindow(width; center=nothing)

Creates a Gaussian window `width` pixels wide, optionally centered at specified position.

# Arguments
- `width`: Window width. Can be a single number or a tuple with width for each dimension.
- `center`: Optional center position. Can be a single number or tuple with center for each dimension.
           If not specified, window is centered at the middle of each dimension.

# Examples
```julia
# Create a symmetric centered Gaussian window
gw = GaussianWindow(15)
window_2d = gw((64, 64))  # centered at (32.5, 32.5)

# Create an asymmetric centered Gaussian window
gw_asym = GaussianWindow((20, 10))
window_asym = gw_asym((64, 32))  # centered at (32.5, 16.5)

# Create an off-center Gaussian window
gw_offset = GaussianWindow(15; center=20)
window_offset = gw_offset((64, 64))  # centered at (20, 20)

# Create an asymmetric off-center window
gw_asym_offset = GaussianWindow((20, 10); center=(30, 15))
window_asym_off = gw_asym_offset((64, 32))  # centered at (30, 15)
```
"""
struct GaussianWindow{T,C} <: Window{T}
    width::T
    center::C

    function GaussianWindow(width::T; center::C=nothing) where {T,C}
        return new{T,C}(width, center)
    end
end

_grange_centered(len, wid) = range(-(len - 1) / wid, (len - 1) / wid; length=len)
_grange_offset(len, wid, center) =
    range((1 - center) / wid, (len - center) / wid; length=len)

_gwidth(w::GaussianWindow{<:Real}, dims) = fill(w.width, length(dims))
_gwidth(w::GaussianWindow{<:Tuple}, dims) =
    if length(w.width) == length(dims)
        w.width
    else
        error("Incompatible dimensions of the window and the array")
    end

_gcenter(w::GaussianWindow{<:Any,Nothing}, dims) = [(d + 1) / 2 for d in dims]  # Default: center of array
_gcenter(w::GaussianWindow{<:Any,<:Real}, dims) = fill(w.center, length(dims))
_gcenter(w::GaussianWindow{<:Any,<:Tuple}, dims) =
    if length(w.center) == length(dims)
        collect(w.center)
    else
        error("Incompatible dimensions of the center and the array")
    end

function (w::GaussianWindow)(dims::Tuple{Vararg{Int}})
    a = zeros(dims) .+ 1
    widths = _gwidth(w, dims)
    centers = _gcenter(w, dims)

    for (id, d) in enumerate(dims)
        ## Create range based on whether center is specified
        if w.center === nothing
            wd = _grange_centered(d, widths[id])
        else
            wd = _grange_offset(d, widths[id], centers[id])
        end

        for (is, s) in enumerate(eachslice(a; dims=id))
            s .*= exp.(-wd[is]^2 / 2)
        end
    end
    return a
end
