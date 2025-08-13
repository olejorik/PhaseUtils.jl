# PhaseUtils.jl

Utilities for phase retrieval workflows. This documentation covers core helpers and shared data structures.

## Quickstart

```julia
using PhaseUtils

# Define a tilt with coefficients [σ, τ₁, τ₂]
t = TiltCentered([0.1, 0.2, -0.3])

# Evaluate on a small Fourier grid
A = materialize(t, (4, 4))
size(A)
```

## Learn more

- Guides
	- [Tilts and Axes](@ref guides/tilts_axes)
- Reference
	- [API](@ref api)

## Index
```@index
```
