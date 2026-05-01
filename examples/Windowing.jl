# ```@meta
# CurrentModule = PhaseUtils
# DocTestSetup = quote
#     using PhaseUtils
# end
# ```

using PhaseUtils
using CairoMakie
CairoMakie.activate!(; type="png")

# # Windowing Functions
#
# Windowing functions are essential tools in signal processing and image analysis for reducing edge effects,
# localizing analysis, and controlling spectral leakage. PhaseUtils provides three window types:
#
# - [`GaussianWindow`](@ref) — smooth Gaussian taper with configurable width and centre position
# - [`HannWindow`](@ref) — raised-cosine taper that reaches exactly 0 at both edges
# - [`TukeyWindow`](@ref) — flat-top taper with controllable ramp fraction `alpha`
#
# All windows are N-dimensional functors: construct a window object, then call it with a size tuple.

# ## Basic Gaussian Windows
#
# Let's start with creating simple centered Gaussian windows of different widths:

# Create windows with different widths
gw_narrow = GaussianWindow(8)
gw_medium = GaussianWindow(15)
gw_wide = GaussianWindow(25)

# Generate 1D windows for visualization
dims_1d = (65,)
window_narrow = gw_narrow(dims_1d)
window_medium = gw_medium(dims_1d)
window_wide = gw_wide(dims_1d)

# In `GaussianWindow`, the `width` parameter is a scaling factor: the 1/√e points fall at
# exactly ±width pixels from the centre of the array.

# Plot 1D windows
fig = Figure(; size=(800, 400))
ax = Axis(
    fig[1, 1];
    title="1D Gaussian Windows with Different Widths",
    xlabel="Pixel Index",
    ylabel="Window Value",
)

lines!(ax, -32:32, window_narrow; label="Width = 8", linewidth=2)
lines!(ax, -32:32, window_medium; label="Width = 15", linewidth=2)
lines!(ax, -32:32, window_wide; label="Width = 25", linewidth=2)

# Reference lines at 1/√e mark the ±width pixel offsets for each window
hlines!(ax, [1 / √ℯ]; color=:red, linestyle=:dot, linewidth=2, alpha=0.8)
vlines!(ax, [-8, 8]; color=:blue, linestyle=:dash, alpha=0.7)
vlines!(ax, [-15, 15]; color=:orange, linestyle=:dash, alpha=0.7)
vlines!(ax, [-25, 25]; color=:green, linestyle=:dash, alpha=0.7)

ax.yticks = ([1 / √ℯ], ["1/√e"])
ax.xticks = ([-25, -15, -8, 8, 15, 25], ["-25", "-15", "-8", "8", "15", "25"])

axislegend(ax; position=:rt)
fig

# ## 2D Gaussian Windows
#
# Gaussian windows are particularly useful in 2D for image processing applications:

# Create 2D windows with symmetric widths
dims_2d = (64, 64)
window_2d_narrow = gw_narrow(dims_2d)
window_2d_medium = gw_medium(dims_2d)
window_2d_wide = gw_wide(dims_2d)

# Visualize 2D windows
fig = Figure(; size=(1000, 300))

ax1 = Axis(fig[1, 1]; title="Narrow (σ=8)", aspect=DataAspect())
heatmap!(ax1, window_2d_narrow; colormap=:viridis)

ax2 = Axis(fig[1, 2]; title="Medium (σ=15)", aspect=DataAspect())
heatmap!(ax2, window_2d_medium; colormap=:viridis)

ax3 = Axis(fig[1, 3]; title="Wide (σ=25)", aspect=DataAspect())
hm = heatmap!(ax3, window_2d_wide; colormap=:viridis)

Colorbar(fig[1, 4], hm; label="Window Value")
fig

# ## Asymmetric Windows
#
# You can create asymmetric windows by specifying different widths for each dimension:

# Create asymmetric windows
gw_horizontal = GaussianWindow((30, 10))  # Wide horizontally, narrow vertically
gw_vertical = GaussianWindow((10, 30))    # Narrow horizontally, wide vertically

window_horizontal = gw_horizontal(dims_2d)
window_vertical = gw_vertical(dims_2d)

fig = Figure(; size=(800, 400))

ax1 = Axis(fig[1, 1]; title="Horizontal Gaussian (30×10)", aspect=DataAspect())
heatmap!(ax1, window_horizontal; colormap=:plasma)

ax2 = Axis(fig[1, 2]; title="Vertical Gaussian (10×30)", aspect=DataAspect())
hm = heatmap!(ax2, window_vertical; colormap=:plasma)

Colorbar(fig[1, 3], hm; label="Window Value")
fig

# ## Off-Center Windows
#
# One of the key features of PhaseUtils' GaussianWindow is the ability to specify custom center positions:

# Create off-center windows
gw_topleft = GaussianWindow(12; center=(16, 16))     # Top-left quadrant
gw_topright = GaussianWindow(12; center=(16, 48))    # Top-right quadrant
gw_bottomleft = GaussianWindow(12; center=(48, 16))  # Bottom-left quadrant
gw_bottomright = GaussianWindow(12; center=(48, 48)) # Bottom-right quadrant

window_tl = gw_topleft(dims_2d)
window_tr = gw_topright(dims_2d)
window_bl = gw_bottomleft(dims_2d)
window_br = gw_bottomright(dims_2d)

fig = Figure(; size=(800, 800))

ax1 = Axis(fig[1, 1]; title="Top-Left (16,16)", aspect=DataAspect())
heatmap!(ax1, window_tl; colormap=:inferno)

ax2 = Axis(fig[1, 2]; title="Top-Right (16,48)", aspect=DataAspect())
heatmap!(ax2, window_tr; colormap=:inferno)

ax3 = Axis(fig[2, 1]; title="Bottom-Left (48,16)", aspect=DataAspect())
heatmap!(ax3, window_bl; colormap=:inferno)

ax4 = Axis(fig[2, 2]; title="Bottom-Right (48,48)", aspect=DataAspect())
hm = heatmap!(ax4, window_br; colormap=:inferno)

Colorbar(fig[:, 3], hm; label="Window Value")
fig

# ## Asymmetric Off-Center Windows
#
# For maximum flexibility, you can combine asymmetric widths with custom center positions:

# Create asymmetric off-center windows
gw_asym_center = GaussianWindow((20, 10); center=(32, 32))    # Centered asymmetric
gw_asym_offset = GaussianWindow((15, 25); center=(20, 45))    # Off-center asymmetric

window_asym_center = gw_asym_center(dims_2d)
window_asym_offset = gw_asym_offset(dims_2d)

fig = Figure(; size=(800, 400))

ax1 = Axis(fig[1, 1]; title="Asymmetric Centered (20×10)", aspect=DataAspect())
heatmap!(ax1, window_asym_center; colormap=:cividis)

ax2 = Axis(fig[1, 2]; title="Asymmetric Off-Center (15×25) at (20,45)", aspect=DataAspect())
hm = heatmap!(ax2, window_asym_offset; colormap=:cividis)

Colorbar(fig[1, 3], hm; label="Window Value")
fig

# ## Cross-Sections and Profiles
#
# Let's examine the profile characteristics of different windows:

# Generate windows for profile analysis
gw_profile = GaussianWindow(15)
gw_offset_profile = GaussianWindow(15; center=(32, 20))

window_profile = gw_profile(dims_2d)
window_offset_profile = gw_offset_profile(dims_2d)

# Extract cross-sections
center_row = window_profile[32, :]
center_col = window_profile[:, 32]
offset_row = window_offset_profile[32, :]
offset_col = window_offset_profile[:, 20]

fig = Figure(; size=(1000, 400))

ax1 = Axis(
    fig[1, 1]; title="Horizontal Profiles", xlabel="Pixel Index", ylabel="Window Value"
)
lines!(ax1, 1:64, center_row; label="Centered", linewidth=2)
lines!(ax1, 1:64, offset_row; label="Off-center", linewidth=2)
axislegend(ax1)

ax2 = Axis(
    fig[1, 2]; title="Vertical Profiles", xlabel="Pixel Index", ylabel="Window Value"
)
lines!(ax2, 1:64, center_col; label="Centered", linewidth=2)
lines!(ax2, 1:64, offset_col; label="Off-center", linewidth=2)
axislegend(ax2)

fig

# ## Hann Window
#
# The [`HannWindow`](@ref) (also called "Hanning") is a raised-cosine taper with no free parameters.
# It reaches exactly 0 at both endpoints and 1 at the centre, making it a good default choice
# when you want guaranteed edge suppression without needing to tune a width parameter.

hw = HannWindow()
window_hann_1d = hw((65,))

fig = Figure(; size=(800, 350))
ax = Axis(
    fig[1, 1];
    title="1D Hann Window",
    xlabel="Pixel Index",
    ylabel="Window Value",
)
lines!(ax, -32:32, window_hann_1d; linewidth=2, color=Makie.wong_colors()[1])
hlines!(ax, [0.0, 1.0]; color=:gray, linestyle=:dot, linewidth=1)
ax.yticks = ([0.0, 0.5, 1.0], ["0", "0.5", "1"])
fig

# The 2D Hann window is the separable outer product of two 1D windows:

window_hann_2d = hw((64, 64))

fig = Figure(; size=(500, 420))
ax = Axis(fig[1, 1]; title="2D Hann Window", aspect=DataAspect())
hm = heatmap!(ax, window_hann_2d; colormap=:viridis)
Colorbar(fig[1, 2], hm; label="Window Value")
fig

# ## Tukey Window
#
# The [`TukeyWindow`](@ref) (tapered cosine) has a single parameter `alpha` that controls
# what fraction of the window is used for the cosine ramps:
#
# - `alpha = 0`: rectangular — no tapering at all
# - `alpha = 1`: identical to a Hann window
# - `0 < alpha < 1`: flat top over the central `1 - alpha` fraction, cosine ramps at each end
#
# This makes it easy to trade off between a wide flat region and strong edge suppression.

alphas = [0.0, 0.25, 0.5, 0.75, 1.0]
N_tukey = 65

fig = Figure(; size=(800, 350))
ax = Axis(
    fig[1, 1];
    title="1D Tukey Windows for Different α",
    xlabel="Pixel Index",
    ylabel="Window Value",
)
for α in alphas
    lines!(ax, -32:32, TukeyWindow(α)((N_tukey,)); label="α = $α", linewidth=2)
end
hlines!(ax, [0.0, 1.0]; color=:gray, linestyle=:dot, linewidth=1)
axislegend(ax; position=:cb)
fig

# The 2D Tukey window for α = 0.5 shows a clear flat-top region in the centre:

window_tukey_2d = TukeyWindow(0.5)((64, 64))

fig = Figure(; size=(500, 420))
ax = Axis(fig[1, 1]; title="2D Tukey Window (α = 0.5)", aspect=DataAspect())
hm = heatmap!(ax, window_tukey_2d; colormap=:viridis)
Colorbar(fig[1, 2], hm; label="Window Value")
fig

# When the two image dimensions require different taper fractions, pass a tuple for `alpha`:

window_tukey_asym = TukeyWindow((0.5, 0.1))((64, 64))

fig = Figure(; size=(500, 420))
ax = Axis(fig[1, 1]; title="2D Tukey Window (α = (0.5, 0.1))", aspect=DataAspect())
hm = heatmap!(ax, window_tukey_asym; colormap=:viridis)
Colorbar(fig[1, 2], hm; label="Window Value")
fig

# ## Window Comparison
#
# Putting all three window types side by side shows the key tradeoffs:
#
# - `GaussianWindow` never reaches 0 — useful when a smooth, infinitely-supported envelope is wanted
# - `HannWindow` guarantees zero edges with no parameters to tune
# - `TukeyWindow(alpha)` lets you preserve a flat signal region while still suppressing edges

N_cmp = 65
pixels = -32:32

fig = Figure(; size=(900, 380))
ax = Axis(
    fig[1, 1];
    title="Window Comparison (N = $N_cmp)",
    xlabel="Pixel Index",
    ylabel="Window Value",
)
lines!(ax, pixels, GaussianWindow(16)((N_cmp,)); label="Gaussian (width=16)", linewidth=2)
lines!(ax, pixels, HannWindow()((N_cmp,));       label="Hann",                linewidth=2)
lines!(ax, pixels, TukeyWindow(0.5)((N_cmp,));   label="Tukey (α=0.5)",      linewidth=2)
lines!(ax, pixels, TukeyWindow(0.25)((N_cmp,));  label="Tukey (α=0.25)",     linewidth=2)
hlines!(ax, [0.0, 1.0]; color=:gray, linestyle=:dot, linewidth=1)
axislegend(ax; position=:cb)
fig

# ## Practical Application: Windowed Fourier Analysis
#
# Applying a window before FFT reduces spectral leakage caused by the implicit periodic extension
# of a finite signal. This example illustrates the resolution–leakage tradeoff: a narrow window
# suppresses background trends but broadens spectral peaks; a wide window preserves peak sharpness
# but retains more background leakage.

using FFTW

## Increase resolution for better spectral analysis
N = 256
x = 1:N
# Add a slowly varying linear background to illustrate edge effects
f1 = 7.3
f2 = 3.7
background = 1 / N .* (x .- N / 2)
test_signal =
    sin.(2π * f1 * x / N) .+ 0.5 * sin.(2π * f2 * x / N) .+ 5background .+ 0.1 * randn(N)

# Create matching windows
window_narrow = GaussianWindow(N / 5)((N,))
window_medium = GaussianWindow(N / 2)((N,))
window_wide = GaussianWindow(N)((N,))

# Apply different windows
windowed_narrow = test_signal .* window_narrow
windowed_medium = test_signal .* window_medium
windowed_wide = test_signal .* window_wide

# Compute FFT spectra, normalised by window energy so peak amplitudes are comparable
fft_original = abs.(fft(test_signal))
fft_narrow = abs.(fft(windowed_narrow))
fft_medium = abs.(fft(windowed_medium))
fft_wide = abs.(fft(windowed_wide))

## Normalise by sqrt(window_energy * N) so a pure sinusoid gives amplitude 0.5 regardless of window
fft_original ./= sqrt(N * N)                              # rectangular implicit window
fft_narrow   ./= sqrt(sum(abs2.(window_narrow)) * N)
fft_medium   ./= sqrt(sum(abs2.(window_medium)) * N)
fft_wide     ./= sqrt(sum(abs2.(window_wide))   * N)

freq = fftfreq(N, 1.0)

f1_freq = f1 / N
f2_freq = f2 / N
## Nearest DFT bin frequencies — spectral peaks land here, not at the true (non-integer) frequencies
f1_bin = round(Int, f1) / N
f2_bin = round(Int, f2) / N

fig = Figure(; size=(1000, 850))

## Time domain
ax1 = Axis(fig[1, 1]; title="Windowed Signals", xlabel="Sample", ylabel="Amplitude")
lines!(ax1, x, test_signal; label="Original", linewidth=2)
lines!(ax1, x, windowed_narrow; label="Narrow window", linewidth=2)
lines!(ax1, x, windowed_medium; label="Medium window", linewidth=2)
lines!(ax1, x, windowed_wide; label="Wide window", linewidth=2)
axislegend(ax1; position=:rt)

## Full frequency domain
ax2 = Axis(fig[2, 1]; title="FFT Spectra", xlabel="Frequency (in 1/N)", ylabel="Magnitude")
lines!(ax2, freq[1:div(N, 2)], fft_original[1:div(N, 2)]; label="Original", linewidth=2)
lines!(ax2, freq[1:div(N, 2)], fft_narrow[1:div(N, 2)]; label="Narrow window", linewidth=2)
lines!(ax2, freq[1:div(N, 2)], fft_medium[1:div(N, 2)]; label="Medium window", linewidth=2)
lines!(ax2, freq[1:div(N, 2)], fft_wide[1:div(N, 2)]; label="Wide window", linewidth=2)
vlines!(ax2, [f2_freq, f1_freq]; linewidth=2, alpha=0.7)
ax2.xticks = ([f2_freq, f1_freq], ["$(f2)", "$(f1)"])
axislegend(ax2; position=:rt)

## Zoom into the peak region
zoom_hw = 2 * f1_freq
ax3 = Axis(
    fig[3, 1];
    title="FFT Spectra — zoom near peaks",
    xlabel="Frequency (in 1/N)",
    ylabel="Magnitude",
    limits=(0, f1_freq + zoom_hw, nothing, nothing),
)
lines!(ax3, freq[1:div(N, 2)], fft_original[1:div(N, 2)]; label="Original", linewidth=2)
lines!(ax3, freq[1:div(N, 2)], fft_narrow[1:div(N, 2)]; label="Narrow window", linewidth=2)
lines!(ax3, freq[1:div(N, 2)], fft_medium[1:div(N, 2)]; label="Medium window", linewidth=2)
lines!(ax3, freq[1:div(N, 2)], fft_wide[1:div(N, 2)]; label="Wide window", linewidth=2)
vlines!(ax3, [f2_freq, f1_freq]; linewidth=2, alpha=0.7)
vlines!(ax3, [f2_bin, f1_bin]; linewidth=1, linestyle=:dash, color=(:gray, 0.5))
ax3.xticks = (
    [f2_bin, f2_freq, f1_bin, f1_freq],
    ["\nk=$(round(Int,f2))/N", "$(f2)", "\nk=$(round(Int,f1))/N", "$(f1)"],
)
axislegend(ax3; position=:rt)

fig


# ## Summary
#
# PhaseUtils provides three complementary window types, all following the same functor interface:
# construct a window object, then call it with a size tuple to get an array.
#
# | Type | Parameters | Reaches 0 at edges? | Flat top? |
# |---|---|---|---|
# | `GaussianWindow(width)` | width (+ optional center) | No | No |
# | `HannWindow()` | none | Yes | No |
# | `TukeyWindow(alpha)` | alpha ∈ \[0,1\] (scalar or tuple) | Yes (for α > 0) | Yes |
#
# All windows are N-dimensional and separable: the N-D window is the outer product of 1D windows
# applied along each dimension. Asymmetric behaviour is available via tuple parameters.
#
# Typical use cases:
# - **Fourier analysis**: apply before FFT to suppress spectral leakage
# - **Image tapering**: suppress edge artefacts before iterative algorithms
# - **Apodisation**: smooth the pupil boundary in phase retrieval
