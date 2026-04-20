# PhaseUtils.jl — Copilot Context

## What this package is

`PhaseUtils` provides small, general-purpose utilities for phase-related computations: wrapping/unwrapping, aperture detection, tilt handling, windowing, and coordinate tools. It is a base-layer dependency used by all other Phase packages and by `Feedback14AMI`.

## Key types and functions

| Symbol | Purpose |
|---|---|
| `phwrap` | Wrap phase to `[-π, π)` |
| `maskedrmse` / `maskedphasermse` | RMS error over aperture mask |
| `find_aperture` | Detect circular aperture from an intensity image |
| `ap2mask` / `mask2ap` | Convert aperture array ↔ boolean mask |
| `circlemask` / `circlemask!` | Generate circular aperture array |
| `Tilt` / `TiltCentered` / `FreeTilt` | Tilt parameterizations: `sigma` (frequency), `tau` (pixel shift) |
| `materialize` / `apply` | Materialize a tilt to an array; apply tilt phase ramp |
| `GaussianWindow` | Apodization window for fringe/PSF processing |
| `ArrayAxes` / `FourierAxes` / `DataAxes` | Coordinate system descriptors |
| `findpiston` / `findtiptilt` | Extract piston/tip-tilt from a phase map |
| `bboxview` | Bounding-box crop of an array |

## Phase unwrapping

Algorithms in `unwrapping.jl` and `algorithms.jl`; includes contour-based and Poisson-solver methods (see `examples/`).

## Relationships

- Used by: `PhaseRetrieval`, `PhaseFromInterferograms`, `PhasePlots`, `Feedback14AMI`
- Depends on: `FFTW`, `FileIO`, `Statistics`

## Source layout

```
src/
    PhaseUtils.jl        ← module entry, exports
    axes_and_tilts.jl    ← Tilt types, ArrayAxes
    aperture_finding.jl  ← find_aperture
    algorithms.jl        ← phase unwrapping algorithms
    windowing.jl         ← GaussianWindow
    utils.jl             ← phwrap, maskedrmse, circlemask, ...
    cropandpad.jl        ← bboxview, padding utilities
```
