## Shared coordinate axes and tilt parameterization for imaging/FT modules

## Dependencies local to this file
using StaticArrays
using LinearAlgebra: dot
using FFTW

## Tilt family
"""
    Tilt

Abstract supertype for affine tilt parameterizations used across imaging
and Fourier-domain modules.

Tilts are represented by a fixed-length coefficient vector `coefs` with the
layout `[σ, τ₁, τ₂, …]`, where `σ` is the constant (piston) term, and
`τ` is the slope vector in the chosen coordinate system. The coordinate
system is provided by an axis policy (see `ArrayAxes` and its subtypes).

Concrete subtypes must store a field `coefs` and are expected to work with
the helpers `sigma`, `tau`, `setsigma!`, `settau!`, `setall!`, `apply`, and
`materialize` defined in this module.
"""
abstract type Tilt end

"""
    sigma(t::Tilt)

Return the constant term `σ` of a tilt.
"""
sigma(t::Tilt) = t.coefs[1]

"""
    tau(t::Tilt)

Return the slope vector `τ` of a tilt.
"""
tau(t::Tilt) = t.coefs[2:end]

"""
    tau(t::Tilt, j)

Return the `j`-th component of the slope vector `τ`.
"""
tau(t::Tilt, j) = t.coefs[1 + j]

"""
    setsigma!(t::Tilt, s)

Set the constant term `σ` in-place.
"""
setsigma!(t::Tilt, s) = (t.coefs[1] = s)

"""
    settau!(t::Tilt, τ)

Set the slope vector `τ` in-place.
"""
settau!(t::Tilt, τ) = (t.coefs[2:end] .= τ)

"""
    setall!(t::Tilt, v)

Set the full coefficient vector `[σ, τ...]` in-place.
"""
setall!(t::Tilt, v) = (t.coefs .= v)

"""
    materialize(t::Tilt, dims)

Build a dense array evaluating the tilt on a grid defined by `dims` using
Fourier-domain coordinates `fftshift(fftfreq(d))` for each dimension `d`.
The result has size `dims`.
"""
materialize(t::Tilt, dims) =
    [sigma(t) + dot(tau(t), x) for x in Iterators.product(fftshift.(fftfreq.(dims))...)]

"""
    materialize(t::Tilt, axes::Vector{<:AbstractVector})

Build a dense array evaluating the tilt over the provided coordinate axes,
one vector per dimension. This is useful with axis policies such as
`FourierAxes()`, `DataAxes()`, or `DataAxesCentered()`.
"""
materialize(t::Tilt, axes::Vector{T} where {T<:AbstractVector}) =
    [sigma(t) + dot(tau(t), x) for x in Iterators.product(axes...)]

"""
    TiltCentered{N} <: Tilt

Concrete tilt storing coefficients in a static vector `coefs::MVector{N,Float64}`
with layout `[σ, τ₁, τ₂, …]`. This is convenient when the tilt is anchored
to a coordinate frame centered at the origin (e.g., frequency coordinates).

Construction: `TiltCentered(coefs::AbstractVector)` converts the input to a
static vector for performance.
"""
struct TiltCentered{N} <: Tilt
    coefs::MVector{N,Float64}
end

TiltCentered(coefs::AbstractVector) = TiltCentered(MVector(coefs...))

"""
    FreeTilt{N} <: Tilt

Concrete tilt with the same storage layout as `TiltCentered`, but without
implied centering semantics. Use this when you want an unconstrained tilt
parameterization in arbitrary coordinates.

Construction: `FreeTilt(coefs::AbstractVector)` converts the input to a
static vector for performance.
"""
struct FreeTilt{N} <: Tilt
    coefs::MVector{N,Float64}
end

FreeTilt(coefs::AbstractVector) = FreeTilt(MVector(coefs...))

apply(t::Tilt, x) = sigma(t) + dot(tau(t), x)
apply(t::Tilt, x, dims) = sum(t.coefs[i + 1] * get(x, i, 1) for i in dims)

## Coordinate axes policies
"""
    ArrayAxes

Abstract supertype for axis policies that produce coordinate vectors for
each dimension of an array. Policies are callable: `axes = policy(dims)`
returns a vector of per-dimension coordinate vectors. These coordinates
define the domain on which tilts are evaluated or materialized.
"""
abstract type ArrayAxes end

"""
    FourierAxes <: ArrayAxes

Axis policy that returns Fourier-domain coordinates for each dimension:
`fftshift(fftfreq(d))`. Useful for representing slopes in normalized
frequency units aligned with FFT conventions.
"""
struct FourierAxes <: ArrayAxes end

"""
    DataAxes <: ArrayAxes

Axis policy that returns 1-based data indices for each dimension: `1:d`.
Useful when slopes are defined with respect to array indices (pixel grid).
"""
struct DataAxes <: ArrayAxes end

"""
    DataAxesCentered <: ArrayAxes

Axis policy that returns centered data coordinates per dimension using
`fftshift(fftfreq(d, d))`, i.e., integer-like coordinates centered at zero.
This is convenient for spatial-domain operations centered on the array.
"""
struct DataAxesCentered <: ArrayAxes end

"""
    ProvidedAxes <: ArrayAxes

Axis policy that wraps explicitly provided coordinate vectors. Use when you
already have precomputed (possibly non-uniform) axes and want to reuse the
same policy interface.

Construction:
    ProvidedAxes(ax1, ax2, ...)
Each axis must be an `AbstractVector`. The policy is callable:
    ProvidedAxes(ax1, ax2)(dims) -> Vector of the stored axes (validates `dims`).
If `dims` do not match the stored axes lengths, an `ArgumentError` is thrown.
"""
struct ProvidedAxes{AX<:Tuple} <: ArrayAxes
    axes::AX
end

ProvidedAxes(axes::AbstractVector...) = ProvidedAxes(tuple(axes...))

function (alg::ProvidedAxes)(dims::NTuple)
    stored_dims = ntuple(i -> length(alg.axes[i]), length(alg.axes))
    stored_dims == dims || throw(
        ArgumentError("ProvidedAxes dims mismatch: got $(dims), expected $(stored_dims)"),
    )
    return [alg.axes...]
end

"""
    (alg::FourierAxes)(dims::NTuple)

Return Fourier-domain coordinate vectors for each dimension in `dims` using
`fftshift(fftfreq(d))`.
"""
(alg::FourierAxes)(dims::NTuple) = [fftshift(fftfreq(d)) for d in dims]
"""
    (alg::DataAxes)(dims::NTuple)

Return 1-based index coordinate vectors `1:d` for each dimension in `dims`.
"""
(alg::DataAxes)(dims::NTuple) = [1:d for d in dims]
"""
    (alg::DataAxesCentered)(dims::NTuple)

Return centered integer-like coordinates for each dimension using
`fftshift(fftfreq(d, d))`.
"""
(alg::DataAxesCentered)(dims::NTuple) = [fftshift(fftfreq(d, d)) for d in dims]

(alg::ArrayAxes)(arr::AbstractArray) = (alg)(size(arr))
