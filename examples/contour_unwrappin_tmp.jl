# ```@meta
# CurrentModule = PhaseUtils
# DocTestSetup = quote
#     using PhaseUtils
#
# using CairoMakie
# CairoMakie.activate!(; type="png")
#     phasedisplay(args...; kwargs...) = heatmap(
#       args...;
#       axis=(aspect=DataAspect(),),
#       colormap=:cyclic_mygbm_30_95_c78_n256,
#       kwargs...,
#       )
#    arraydisplay(args...; kwargs...) = heatmap(args...; axis=(aspect=DataAspect(),), kwargs...)
# # end
# ```
using PhaseUtils
using CairoMakie
CairoMakie.activate!(; type="png")
phasedisplay(args...; kwargs...) = heatmap(
    args...;
    axis=(aspect=DataAspect(),),
    colormap=:cyclic_mygbm_30_95_c78_n256,
    kwargs...,
)
arraydisplay(args...; kwargs...) = heatmap(args...; axis=(aspect=DataAspect(),), kwargs...)

# # Temporary file, for test/debug purposes

# ##  Example with circular aperture (for tests)
s1, s2, m = 140, 100, 25
ap = zeros(s1, s2)
y = range(-1.1 * s1 / s2, 1.1 * s1 / s2, s1)
x = range(-1.1, 1.1, s2)
ap[[x .^ 2 + y .^ 2 .<= 1 for y in y, x in x]] .= 1
phaseGT = [-x^3 + 3x .^ 2 + y .^ 2 - 10y for y in y, x in x]
phaseGT .-= sum(phaseGT .* ap) / sum(ap) #extract mean
phase = phwrap(phaseGT) .* ap
phasedisplay(phase)

sol = unwrap_LS(phase, ap; restore_piston=false)
sol .-= sum(sol .* ap) / sum(ap) #extract mean
fig, ax, hm = arraydisplay((sol .- phaseGT) .* ap2mask(ap))
ax.title = "Phase restoration error is $(maskedrmse(sol, phaseGT, ap)) rms"
Colorbar(fig[1, 2])
fig
