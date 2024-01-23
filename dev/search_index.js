var documenterSearchIndex = {"docs":
[{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"EditURL = \"../../../examples/Unwrapping.jl\"","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"CurrentModule = PhaseUtils\nDocTestSetup = quote\n    using PhaseUtils\nend","category":"page"},{"location":"examples/Unwrapping/#Phase-Unwrapping","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"","category":"section"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"We will use the calculated response function for the phase unwrapping demonstration.","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"using PhaseUtils\nusing CairoMakie\nCairoMakie.activate!(; type=\"png\")\nmask = circlemask((300, 300), 150, 150, 145)\nact = circlemask((300, 300), 50, 153.5, 10)\nresp = membrane_sor(act, mask);\nnothing #hide","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"First, we define a couple of plotting functions to save on typing","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"phasedisplay(args...; kwargs...) = heatmap(\n    args...;\n    axis=(aspect=DataAspect(),),\n    colormap=:cyclic_mygbm_30_95_c78_n256,\n    kwargs...,\n)\narraydisplay(args...; kwargs...) = heatmap(args...; axis=(aspect=DataAspect(),), kwargs...)\n\nresp ./= 10\narraydisplay(resp)\n\nresp_wr = phwrap(resp)\nphasedisplay(resp_wr)","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"This noiseless, periodic, and defined on the whole domain phase is easy to unwrap by line-by-line integration of wrapped phase differences using the Itoh algorithm, see itoh.","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"ph_itoh = reshape(itoh(resp_wr[:]), size(resp_wr))\narraydisplay(ph_itoh)","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"Phases occurring in real life are usually a little bit more complicated. First of all, they are often defined only in some region Ω, and do not zero out on the boundary of the region Ω`. This is simulated below by adding some linear tilt to the phase.","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"tilt = [(i[1] * 0.000135 + i[2] * 0.0011 + 10) for i in eachindex(IndexCartesian(), mask)]\nresp_wr = phwrap((resp .+ tilt)) .* mask\nphasedisplay(resp_wr .* ap2mask(mask))","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"(here we have used ap2mask function which removes the pixels outside the mask for better visibility).","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"Here is what happens if we subtract the (exactly known) tilt from our new wrapped phase","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"phasedisplay((resp_wr .- tilt) .* ap2mask(mask))","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"This phase should not be restorable by the Itoh algorithm because the jumps on the boundary are greater than π, and indeed we see the errors in unwrapping","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"ph_itoh = reshape(itoh(resp_wr[:]), size(resp_wr))\narraydisplay((ph_itoh .- tilt) .* mask)","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"but wrapped back it looks OK, of course (what we see is just unwrapping errors, assigning wrong k to an unwrapped value hatphi = psi + 2πk k  mathbbZ):","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"phasedisplay(phwrap(ph_itoh .- tilt) .* ap2mask(mask))","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"We can unwrap it with the least-squares algorithm","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"phi_LS = unwrap_LS(resp_wr, mask)\narraydisplay((phi_LS .- tilt) .* ap2mask(mask))","category":"page"},{"location":"examples/Unwrapping/#Example-with-wedge","page":"Phase Unwrapping","title":"Example with wedge","text":"","category":"section"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"wedge = mask .* linearphase((300, 300), 150, 150, 0.0, 0.2)\nfig, ax, hm = arraydisplay(wedge)\nColorbar(fig[1, 2], hm)\nax.title = \"Original wedge\"\nfig\n\nwedge_wr = phwrap(wedge)\nphasedisplay(wedge_wr)","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"It's impossible to unwrap it with the Itoh algorithm.","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"wedge_itoh = reshape(itoh(wedge_wr[:]), size(wedge_wr))\nfig, ax, hm = arraydisplay(wedge_itoh)\nColorbar(fig[1, 2], hm)\nax.title = \"RMS error is $(maskedrmse(wedge_itoh, wedge, mask))\"\nfig","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"Nor it's possible to do this using Itoh's algorithm along the rows.","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"wedge_itoh = reshape(itoh(wedge_wr'[:]), size(wedge_wr))'\nfig, ax, hm = arraydisplay(wedge_itoh)\nColorbar(fig[1, 2], hm)\nax.title = \"RMS error is $(maskedrmse(wedge_itoh, wedge, mask))\"\nfig","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"But we can unwrap it with the least-squares algorithm","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"wedge_LS = unwrap_LS(wedge_wr, mask)\nfig, ax, hm = arraydisplay(wedge_LS)\nColorbar(fig[1, 2], hm)\nfig","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"Compare with the original wedge:","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"fig, ax, hm = arraydisplay((wedge_LS .- wedge))\nColorbar(fig[1, 2], hm)\nax.title = \"RMS error is $(maskedrmse(wedge_LS, wedge, mask))\"\nfig","category":"page"},{"location":"examples/Unwrapping/#Unwrapping-of-the-phase-with-residues","page":"Phase Unwrapping","title":"Unwrapping of the phase with residues","text":"","category":"section"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"Let us add some bad pixels to the wrapped wedge","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"resblock!(arr, r, c) = (arr[c:(c + 1), r:(r + 1)] .= [0 π/2; -π/2 π])\nresblock!(wedge_wr, 100, 110)\nresblock!(wedge_wr, 50, 150)\nresblock!(wedge_wr, 130, 130)\n\nfig, ax, hm = phasedisplay(wedge_wr)\nax.title = \"RMS error with the original phase is $(maskedrmse(wedge_wr, phwrap(wedge), mask))\"\nfig","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"It's still impossible to unwrap it with the Itoh algorithm, but now the error is of a different type","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"wedge_itoh = reshape(itoh(wedge_wr[:]), size(wedge_wr))\nfig, ax, hm = arraydisplay(wedge_itoh)\nax.title = \"RMS error is $(maskedrmse(wedge_itoh, wedge, mask))\"\nColorbar(fig[1, 2], hm)\nfig","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"Nor it's possible to do this using Itoh's algorithm along the rows.","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"wedge_itoh = reshape(itoh(wedge_wr'[:]), size(wedge_wr))'\nfig, ax, hm = arraydisplay(wedge_itoh)\nax.title = \"RMS error is $(maskedrmse(wedge_itoh, wedge, mask))\"\nColorbar(fig[1, 2], hm)\nfig","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"But we can unwrap it with the least-squares algorithm","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"wedge_LS = unwrap_LS(wedge_wr, mask)\nfig, ax, hm = arraydisplay(wedge_LS)\nColorbar(fig[1, 2], hm)\nfig","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"Compare with the original wedge:","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"fig, ax, hm = arraydisplay(ap2mask(mask) .* (wedge_LS .- wedge))\nax.title = \"RMS error is $(maskedrmse(wedge_LS, wedge, mask))\"\nColorbar(fig[1, 2], hm)\nfig","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"","category":"page"},{"location":"examples/Unwrapping/","page":"Phase Unwrapping","title":"Phase Unwrapping","text":"This page was generated using Literate.jl.","category":"page"},{"location":"examples/contour_unwrappin_tmp/","page":"Temporary file, for test/debug purposes","title":"Temporary file, for test/debug purposes","text":"EditURL = \"../../../examples/contour_unwrappin_tmp.jl\"","category":"page"},{"location":"examples/contour_unwrappin_tmp/","page":"Temporary file, for test/debug purposes","title":"Temporary file, for test/debug purposes","text":"CurrentModule = PhaseUtils\nDocTestSetup = quote\n    using PhaseUtils\n\nusing CairoMakie\nCairoMakie.activate!(; type=\"png\")\n    phasedisplay(args...; kwargs...) = heatmap(\n      args...;\n      axis=(aspect=DataAspect(),),\n      colormap=:cyclic_mygbm_30_95_c78_n256,\n      kwargs...,\n      )\n   arraydisplay(args...; kwargs...) = heatmap(args...; axis=(aspect=DataAspect(),), kwargs...)\n# end","category":"page"},{"location":"examples/contour_unwrappin_tmp/","page":"Temporary file, for test/debug purposes","title":"Temporary file, for test/debug purposes","text":"using PhaseUtils\nusing CairoMakie\nCairoMakie.activate!(; type=\"png\")\nphasedisplay(args...; kwargs...) = heatmap(\n    args...;\n    axis=(aspect=DataAspect(),),\n    colormap=:cyclic_mygbm_30_95_c78_n256,\n    kwargs...,\n)\narraydisplay(args...; kwargs...) = heatmap(args...; axis=(aspect=DataAspect(),), kwargs...)","category":"page"},{"location":"examples/contour_unwrappin_tmp/#Temporary-file,-for-test/debug-purposes","page":"Temporary file, for test/debug purposes","title":"Temporary file, for test/debug purposes","text":"","category":"section"},{"location":"examples/contour_unwrappin_tmp/#Example-with-circular-aperture-(for-tests)","page":"Temporary file, for test/debug purposes","title":"Example with circular aperture (for tests)","text":"","category":"section"},{"location":"examples/contour_unwrappin_tmp/","page":"Temporary file, for test/debug purposes","title":"Temporary file, for test/debug purposes","text":"s1, s2, m = 140, 100, 25\nap = zeros(s1, s2)\ny = range(-1.1 * s1 / s2, 1.1 * s1 / s2, s1)\nx = range(-1.1, 1.1, s2)\nap[[x .^ 2 + y .^ 2 .<= 1 for y in y, x in x]] .= 1\nphaseGT = [-x^3 + 3x .^ 2 + y .^ 2 - 10y for y in y, x in x]\nphaseGT .-= sum(phaseGT .* ap) / sum(ap) #extract mean\nphase = phwrap(phaseGT) .* ap\nphasedisplay(phase)\n\nsol = unwrap_LS(phase, ap; restore_piston=false)\nsol .-= sum(sol .* ap) / sum(ap) #extract mean\nfig, ax, hm = arraydisplay((sol .- phaseGT) .* ap2mask(ap))\nax.title = \"Phase restoration error is $(maskedrmse(sol, phaseGT, ap)) rms\"\nColorbar(fig[1, 2])\nfig","category":"page"},{"location":"examples/contour_unwrappin_tmp/","page":"Temporary file, for test/debug purposes","title":"Temporary file, for test/debug purposes","text":"","category":"page"},{"location":"examples/contour_unwrappin_tmp/","page":"Temporary file, for test/debug purposes","title":"Temporary file, for test/debug purposes","text":"This page was generated using Literate.jl.","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"EditURL = \"../../../examples/Poisson.jl\"","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"CurrentModule = PhaseUtils\nDocTestSetup = quote\n    using PhaseUtils\nend","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"using PhaseUtils\nusing CairoMakie\nCairoMakie.activate!(; type=\"png\")","category":"page"},{"location":"examples/Poisson/#Solving-Poisson-equation","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"","category":"section"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"Use membrane_sor(f, a) to calculate the solution of the Poisson equation with zero boundary conditions","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"beginaligned\nDelta u(x) =  f(x) quad x in Ω\nu(x) =   0 quad  x in partial Ω\nendaligned","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"Region Ω is defined by the boolean array a, which should have the same dimensions as the source array f.","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"Let us simulate a response of a membrane deformable mirror using the Poisson equation. Define a circular aperture of 145-pixel radius and a circular actuator  of 10-pixel radius inside it:","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"mask = circlemask((300, 300), 150, 150, 145)\nact = circlemask((300, 300), 50, 153.5, 10);\nnothing #hide","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"Plot actuator inside the aperture:","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"heatmap(mask .+ act; axis=(aspect=DataAspect(),))","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"Calculate the response and show the results","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"resp = membrane_sor(act, mask)\ncontourf(resp; axis=(aspect=DataAspect(),))","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"Check that the Laplacian of the calculated response is proportional to the actuator shape:","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"act_restored = PhaseUtils._calculate_Laplacian(resp)\nheatmap(act_restored .* mask; axis=(aspect=DataAspect(),))","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"To set the non-zero boundary conditions, use the mutating version of the solver:","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"maskrect = zeros(Bool, 300, 300)\nmaskrect[2:(end - 1), 2:(end - 1)] .= true\nu = zeros(size(mask))\nu[1, :] .= 1\nu[end, 100:200] .= -1\nmembrane_sor!(u, zeros(size(u)), maskrect)\nfig, ax, hm = heatmap(u; axis=(aspect=DataAspect(),))\ncontour!(u; labels=true, levels=-1:0.1:1, labelsize=15, color=:black)\nfig","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"Next example sets the boundary values for a circular aperture","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"edge_out, _ = find_cw_border(mask; outside=true)\nu = zeros(size(mask))\nu[edge_out] .= sin.(range(0, 4 * 2π, length(edge_out) + 1)[1:(end - 1)])\nheatmap(u; axis=(aspect=DataAspect(),))\n\nmembrane_sor!(u, zeros(size(u)), mask; maxits=1000)\nfig, ax, hm = heatmap(u .* ap2mask(mask); axis=(aspect=DataAspect(),))\ncontour!(u; labels=true, levels=-0.9:0.1:0.9, labelsize=15, color=:white)\nfig","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"","category":"page"},{"location":"examples/Poisson/","page":"Solving Poisson equation","title":"Solving Poisson equation","text":"This page was generated using Literate.jl.","category":"page"},{"location":"about/#About-the-package","page":"About","title":"About the package","text":"","category":"section"},{"location":"about/","page":"About","title":"About","text":"This package contains small utilities used often in Phase-Retrieval-related context.","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"EditURL = \"../../../examples/contours_detection.jl\"","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"CurrentModule = PhaseUtils\nDocTestSetup = quote\n    using PhaseUtils\nend","category":"page"},{"location":"examples/contours_detection/#Test-of-the-contour-detection","page":"Test of the contour detection","title":"Test of the contour detection","text":"","category":"section"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"I use simple contour detection for phase unwrapping. The goal is to find the points where the finite difference is not defined.","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"Begin with a circular aperture.","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"using PhaseUtils\nusing CairoMakie\nCairoMakie.activate!(; type=\"png\")\narraydisplay(args...; kwargs...) = heatmap(args...; axis=(aspect=DataAspect(),), kwargs...)\n\nap = circlemask((30, 20), 10, 10, 8)\nfig, ax, hm = arraydisplay(ap; colormap=:reds, alpha=0.25)\nedg = PhaseUtils._find_set_edges(ap)\nfor k in keys(edg)\n    scatter!(\n        map(x -> Tuple(x) .+ Tuple(PhaseUtils.step[k]) .* 0.2, edg[k]); label=\"$(String(k))\"\n    )\nend\nleg = axislegend(ax)\nax.title = \"Detected edge pixels for inner circle\"\nfig","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"We can now sort all these pixels in a clockwise manner","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"cwborder, _ = PhaseUtils._find_cw_border_alloc(ap)\nscatter!(map(Tuple, cwborder); marker='o', color=:black, label=\"Contour\")\ndelete!(leg)\naxislegend(ax)\nfig","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"Now check the same for the inverse mask","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"notap = .!ap\nfig, ax, hm = arraydisplay(notap; colormap=:reds, alpha=0.25)\nedg = PhaseUtils._find_set_edges(notap)\nfor k in keys(edg)\n    scatter!(\n        map(x -> Tuple(x) .+ Tuple(PhaseUtils.step[k]) .* 0.2, edg[k]); label=\"$(String(k))\"\n    )\nend\nleg = axislegend(ax)\nax.title = \"Detected edge pixels for outer circle\"\nfig","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"We can now sort all these pixels in a clockwise manner","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"cwborder, _ = PhaseUtils._find_cw_border_alloc(ap; outside=true)\nscatter!(map(Tuple, cwborder); marker='o', color=:black, label=\"Contour\")\ndelete!(leg)\naxislegend(ax)\nfig","category":"page"},{"location":"examples/contours_detection/#More-difficult-shapes","page":"Test of the contour detection","title":"More difficult shapes","text":"","category":"section"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"ap = circlemask((40, 22), 10, 11, 8) .|| circlemask((40, 22), 26, 11, 9)\nfig, ax, hm = arraydisplay(ap; colormap=:reds, alpha=0.25)\nedg = PhaseUtils._find_set_edges(ap)\nfor k in keys(edg)\n    scatter!(\n        map(x -> Tuple(x) .+ Tuple(PhaseUtils.step[k]) .* 0.2, edg[k]); label=\"$(String(k))\"\n    )\nend\nleg = axislegend(ax)\nax.title = \"Detected edge pixels for inner edge\"\nfig","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"We can now sort all this pixels in a clockwise manner","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"cwborder, _ = PhaseUtils._find_cw_border_alloc(ap)\nscatter!(map(Tuple, cwborder); marker='o', color=:black, label=\"Contour\")\ndelete!(leg)\naxislegend(ax)\nfig","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"Now check the same for the inverse mask","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"notap = .!ap\nfig, ax, hm = arraydisplay(notap; colormap=:reds, alpha=0.25)\nedg = PhaseUtils._find_set_edges(notap)\nfor k in keys(edg)\n    scatter!(\n        map(x -> Tuple(x) .+ Tuple(PhaseUtils.step[k]) .* 0.2, edg[k]); label=\"$(String(k))\"\n    )\nend\nleg = axislegend(ax)\nax.title = \"Detected edge pixels for outer edge\"\nfig","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"We can now sort all these pixels in a clockwise manner","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"cwborder, _ = PhaseUtils._find_cw_border_alloc(ap; outside=true)\nscatter!(map(Tuple, cwborder); marker='o', color=:black, label=\"Contour\")\ndelete!(leg)\naxislegend(ax)\nfig","category":"page"},{"location":"examples/contours_detection/#Conclusions","page":"Test of the contour detection","title":"Conclusions","text":"","category":"section"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"It works as expected.","category":"page"},{"location":"examples/contours_detection/#Dual-and-half-dual-grids","page":"Test of the contour detection","title":"Dual and half-dual grids","text":"","category":"section"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"Fro Talmi-Ribak's method, I also need to find points on the dual grid that fall within the area and on its boundary (with at least two adjacent pixels). Here is a straightforward approach.","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"function dual_region(ap)\n    dualap = zeros(Bool, size(ap) .+ 1) # index i,j corresponds to the point i-1/2,j-1/2\n    for ind in eachindex(IndexCartesian(), ap)\n        if ap[ind]\n            if ap[ind + CartesianIndex(1, 0)] #two horisontal neighbours (i,j) and (i+1,j)\n                dualap[ind + CartesianIndex(1, 0)] = 1 # (i+1/2,j-1/2)\n                dualap[ind + CartesianIndex(1, 1)] = 1 # (i+1/2,j+1/2)\n            end\n            if ap[ind + CartesianIndex(0, 1)] #two vertical neighbours (i,j) and (i,j+1)\n                dualap[ind + CartesianIndex(0, 1)] = 1 # (i-1/2,j+1/2)\n                dualap[ind + CartesianIndex(1, 1)] = 1 # (i+1/2,j+1/2)\n            end\n        end\n    end\n    return dualap\nend\n\nfunction dualcoordmin(ind)\n    return ind.I .- 0.5\nend\nfunction dualcoordplus(ind)\n    return ind.I .+ 0.5\nend","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"Plot the ap and its dual","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"fig, ax, hm = arraydisplay(ap; colormap=:reds, alpha=0.25);\ndualap = dual_region(ap)\ncwborder, _ = PhaseUtils._find_cw_border_alloc(dualap)\nscatter!(dualcoordmin.(findall(dualap)); marker='+', color=:blue, label=\"dual region\")\nscatter!(map(dualcoordmin, cwborder); color=:blue, marker='□', label=\"dual contour\")\nleg = axislegend(ax)\nax.title = \"Detected dual region and its border\"\nfig","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"","category":"page"},{"location":"examples/contours_detection/","page":"Test of the contour detection","title":"Test of the contour detection","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#PhaseUtils.jl","page":"Home","title":"PhaseUtils.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for PhaseUtils.jl","category":"page"},{"location":"#Types-and-Functions","page":"Home","title":"Types and Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [PhaseUtils]\nPrivate = false","category":"page"},{"location":"#PhaseUtils.ap2mask-Tuple{Any}","page":"Home","title":"PhaseUtils.ap2mask","text":"ap2mask(ap)\n\nConverts  ap array to mask, so that every zrro is mapped to NaN and non-zero elements to 1.\n\nSee also mask2ap.\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.bboxview","page":"Home","title":"PhaseUtils.bboxview","text":"bboxview(arr, mask, pad = 0)\n\nMake a box corresponding to the mask not-NaN elements surrounded by pad.\n\n\n\n\n\n","category":"function"},{"location":"#PhaseUtils.circlemask!-Tuple{Matrix{<:Number}, Any, Any, Any}","page":"Home","title":"PhaseUtils.circlemask!","text":"circlemask!(a::Matrix{<:Number}, [cx, cy,] r)\n\nSet to zero values of matrix a outside the circle with center (cx,cy) and radius r.\n\nIf omitted, (cx,cy) is the center of the matrix.\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.circlemask-Tuple{Tuple{Int64, Int64}, Vararg{Any}}","page":"Home","title":"PhaseUtils.circlemask","text":"circlemask(dims::NTuple{2, Int}, cx, cy, r)\n\nCreate a boolean matrix of size dims with true only inside the circle with center (cx,cy) and radius r.\n\nIf omitted, (cx,cy) is the center of the domain: (dims+1)/2.\n\nExamples\n\njulia> circlemask((6,8), 2.5,3, 1.5)\n6×8 Matrix{Bool}:\n 0  0  1  0  0  0  0  0\n 0  1  1  1  0  0  0  0\n 0  1  1  1  0  0  0  0\n 0  0  1  0  0  0  0  0\n 0  0  0  0  0  0  0  0\n 0  0  0  0  0  0  0  0\n\njulia> circlemask((10,10), 3)\n10×10 Matrix{Bool}:\n 0  0  0  0  0  0  0  0  0  0\n 0  0  0  0  0  0  0  0  0  0\n 0  0  0  1  1  1  1  0  0  0\n 0  0  1  1  1  1  1  1  0  0\n 0  0  1  1  1  1  1  1  0  0\n 0  0  1  1  1  1  1  1  0  0\n 0  0  1  1  1  1  1  1  0  0\n 0  0  0  1  1  1  1  0  0  0\n 0  0  0  0  0  0  0  0  0  0\n 0  0  0  0  0  0  0  0  0  0\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.crop-Tuple{AbstractArray, Any, Any}","page":"Home","title":"PhaseUtils.crop","text":"crop(a::AbstractArray, cropsize[, cropcenter])\ncrop(cropsize[, cropcenter]) creates a cropping operator.\n\nCrop the array to cropsize around cropcenter. If cropcenter is omitted, use central element of the array.\n\ncropsize and cropcenter can be either tuples/vectors of integers or just single integer, which represent a tuple of identical numbers.\n\nExamples\n\njulia> a = reshape(1:15, 3, 5); crop(a, 1, (1, 3))\n1×1 Matrix{Int64}:\n 7\n\njulia> crop(a, (2, 2), (1, 3))\n1×2 Matrix{Int64}:\n 4  7\n\njulia> crop(3)(a)\n3×3 Matrix{Int64}:\n 4  7  10\n 5  8  11\n 6  9  12\n\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.integrate_2dgrad","page":"Home","title":"PhaseUtils.integrate_2dgrad","text":"integrate_2dgrad(gx, gy[, gradmethod=default_grad_method(gx)])\n\nTBW\n\n\n\n\n\n","category":"function"},{"location":"#PhaseUtils.integrate_periodic_grad-Tuple{Any}","page":"Home","title":"PhaseUtils.integrate_periodic_grad","text":"integrate_periodic_grad(g)\n\nIntegrate 1D gradient with periodic boundary conditions\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.itoh-Tuple{Any}","page":"Home","title":"PhaseUtils.itoh","text":"itoh(phi)\n\nItoh's algorithm for 1D phase unwrapping.\n\nExample\n\njulia> ph = 1:10\n1:10\n\njulia> ph_w = phwrap(ph)\n10-element Vector{Float64}:\n  1.0\n  2.0\n  3.0\n -2.2831853071795867\n -1.2831853071795865\n -0.28318530717958645\n  0.7168146928204135\n  1.7168146928204135\n  2.7168146928204133\n -2.566370614359173\n\njulia> itoh(ph_w)\n10-element Vector{Float64}:\n  1.0\n  2.0\n  3.0\n  4.0\n  5.0\n  6.0\n  7.0\n  8.0\n  9.0\n 10.0\n\n\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.linearphase-Tuple{Type{<:Number}, Tuple{Int64, Int64}, Vararg{Any, 4}}","page":"Home","title":"PhaseUtils.linearphase","text":"linearphase(T=Float64, dims::NTuple{2,Int},  cx, cy, kx, ky)\n\nCreate a matrix with linear phase with slopes (kx, ky) and having zero at point (cx, cy).\n\nExample\n\njulia> linearphase((5,7), 3, 3, 1.5, 1.5)\n5×7 Matrix{Float64}:\n -6.0  -4.5  -3.0  -1.5  0.0  1.5  3.0\n -4.5  -3.0  -1.5   0.0  1.5  3.0  4.5\n -3.0  -1.5   0.0   1.5  3.0  4.5  6.0\n -1.5   0.0   1.5   3.0  4.5  6.0  7.5\n  0.0   1.5   3.0   4.5  6.0  7.5  9.0\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.mask2ap-Tuple{Any}","page":"Home","title":"PhaseUtils.mask2ap","text":"mask2ap(mask)\n\nConverts  mask array where points outside the aperture are defines as NaN to a Float array, so that NaN -> 0, not NaN -> 1.\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.membrane_sor!-Tuple{Any, Array, Matrix{Bool}}","page":"Home","title":"PhaseUtils.membrane_sor!","text":"membrane_sor!(u, f::Array{Real}, mask_internal::Array{Bool})\n\nSolves Poisson equation Δ u = f, with boundary conditions set by u(x) for x ∉ mask_internal using Successive Overrelaxation method (see Numerical Recipes (www.cambridge.org/9780521880688), ch.20.5.1).\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.membrane_sor-Tuple{Array, Matrix{Bool}}","page":"Home","title":"PhaseUtils.membrane_sor","text":"membrane_sor(f::Array{Real}, mask_internal::Array{Bool})\n\nSolves Poisson equation Δ u = f, with boundary conditions u(x) = 0 <=> x ∉ mask_internal using Successive Overrelaxation method (see Numerical Recipes (www.cambridge.org/9780521880688), ch.20.5.1).\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.phwrap-Tuple{Number}","page":"Home","title":"PhaseUtils.phwrap","text":"phwrap(ϕ)\n\nWrap phase ϕ: ψ = phwrap(ϕ) ⇔ ψ = ϕ + 2πk, k∈Z, -π <ψ ≤π. If called on an array, works pointwise.\n\n\n\n\n\n","category":"method"},{"location":"#PhaseUtils.unwrap_LS-Tuple{Any, Any}","page":"Home","title":"PhaseUtils.unwrap_LS","text":"unwrap_LS(phase, aperture; restore_piston=true)\n\nUnwrap 2D phase defined inside aperture using the Least-Squares decomposition of the wrapped gradient of the wprapped phase in the rotor-free and solenodial field and integration of the rotor-free part.\n\n\n\n\n\n","category":"method"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
