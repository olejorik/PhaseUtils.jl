"""
    Windowing functionality for Phase project

This module provides windowing functions for image processing and signal analysis.
"""

abstract type Window end

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
struct GaussianWindow{T,C} <: Window
    width::T
    center::C

    function GaussianWindow(width::T; center::C=nothing) where {T,C}
        return new{T,C}(width, center)
    end
end

# Single canonical range: runs from (1 - center)/wid to (len - center)/wid.
# At the default center (len+1)/2 this is symmetric around zero with half-range (len-1)/(2*wid).
_grange(len, wid, center) = range((1 - center) / wid, (len - center) / wid; length=len)

_gwidth(w::GaussianWindow{<:Real}, dims) = fill(w.width, length(dims))
_gwidth(w::GaussianWindow{<:Tuple}, dims) =
    if length(w.width) == length(dims)
        collect(w.width)
    else
        error("Incompatible dimensions of the window and the array")
    end

_gcenter(w::GaussianWindow{<:Any,Nothing}, dims) = [(d + 1) / 2 for d in dims]
_gcenter(w::GaussianWindow{<:Any,<:Real}, dims) = fill(w.center, length(dims))
_gcenter(w::GaussianWindow{<:Any,<:Tuple}, dims) =
    if length(w.center) == length(dims)
        collect(w.center)
    else
        error("Incompatible dimensions of the center and the array")
    end

function (w::GaussianWindow)(dims::Tuple{Vararg{Int}})
    a = ones(dims)
    widths = _gwidth(w, dims)
    centers = _gcenter(w, dims)

    for (id, d) in enumerate(dims)
        wd = _grange(d, widths[id], centers[id])

        for (is, s) in enumerate(eachslice(a; dims=id))
            s .*= exp.(-wd[is]^2 / 2)
        end
    end
    return a
end

# ---------------------------------------------------------------------------
# HannWindow
# ---------------------------------------------------------------------------

"""
    HannWindow()

Hann (raised-cosine) window, also commonly called "Hanning".

The 1D window of length `N` is defined as
```
w[i] = 0.5 * (1 - cos(2π(i-1)/(N-1))),  i = 1 … N
```
so it is exactly 0 at both endpoints and 1 at the centre.

The N-D window is the separable outer product of 1D windows along each dimension.

# Example
```julia
hw = HannWindow()
w1d = hw((64,))     # 64-point 1D Hann window
w2d = hw((64, 64))  # 64×64 separable 2D Hann window
```
"""
struct HannWindow <: Window end

function _hann1d(n::Int)
    return [0.5 * (1 - cospi(2 * (i - 1) / (n - 1))) for i in 1:n]
end

function (w::HannWindow)(dims::Tuple{Vararg{Int}})
    a = ones(dims)
    for (id, d) in enumerate(dims)
        wd = _hann1d(d)
        for (is, s) in enumerate(eachslice(a; dims=id))
            s .*= wd[is]
        end
    end
    return a
end

# ---------------------------------------------------------------------------
# TukeyWindow
# ---------------------------------------------------------------------------

"""
    TukeyWindow(alpha)

Tukey (tapered cosine) window with taper fraction `alpha`.

- `alpha = 0`: rectangular window (no taper).
- `alpha = 1`: equivalent to a Hann window.
- `0 < alpha < 1`: flat top of fractional width `1 - alpha`, with cosine
  ramps of fractional width `alpha/2` at each end.

`alpha` can be a single number (applied to all dimensions) or a tuple to
set a different taper fraction per dimension.

The 1D window of length `N` with taper fraction `α` is:
```
w[i] = 0.5 * (1 - cos(π(i-1) / (αN/2 - 1))),  1 ≤ i ≤ αN/2
w[i] = 1,                                        αN/2 < i < N(1 - α/2)
w[i] = 0.5 * (1 - cos(π(N-i)   / (αN/2 - 1))), N(1 - α/2) ≤ i ≤ N
```

The N-D window is the separable outer product of 1D windows along each dimension.

# Example
```julia
tw = TukeyWindow(0.5)
w1d = tw((64,))       # 64-point 1D Tukey window, 50 % taper
w2d = tw((64, 64))    # 64×64 separable 2D Tukey window

tw_asym = TukeyWindow((0.5, 0.1))
w2d_asym = tw_asym((64, 64))   # different taper per dimension
```
"""
struct TukeyWindow{T} <: Window
    alpha::T

    function TukeyWindow(alpha::T) where {T<:Union{Real,Tuple}}
        return new{T}(alpha)
    end
end

_talpha(w::TukeyWindow{<:Real}, dims) = fill(w.alpha, length(dims))
_talpha(w::TukeyWindow{<:Tuple}, dims) =
    if length(w.alpha) == length(dims)
        collect(w.alpha)
    else
        error("Incompatible dimensions of TukeyWindow alpha and the array")
    end

function _tukey1d(n::Int, α::Real)
    α = clamp(α, 0, 1)
    w = ones(n)
    # Standard Tukey definition: half-ramp spans α*(N-1)/2 samples.
    # At α=1 this reduces to the Hann formula exactly (for any N).
    half_ramp = α * (n - 1) / 2
    ramp = floor(Int, half_ramp) + 1      # number of left ramp samples
    ramp < 2 && return w                  # α ≈ 0: rectangular
    for i in 1:ramp
        v = 0.5 * (1 - cospi((i - 1) / half_ramp))
        w[i] = v
        w[n - i + 1] = v
    end
    return w
end

function (w::TukeyWindow)(dims::Tuple{Vararg{Int}})
    a = ones(dims)
    alphas = _talpha(w, dims)
    for (id, d) in enumerate(dims)
        wd = _tukey1d(d, alphas[id])
        for (is, s) in enumerate(eachslice(a; dims=id))
            s .*= wd[is]
        end
    end
    return a
end
