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
# localizing analysis, and controlling spectral leakage. PhaseUtils provides flexible Gaussian windowing
# functionality through the [`GaussianWindow`](@ref) type.
#
# The `GaussianWindow` allows you to create smooth, bell-shaped windows with configurable width and center position,
# making it ideal for applications in phase analysis, interferometry, and Fourier domain processing.

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

# Note: In our GaussianWindow implementation, the `width` parameter controls how the Gaussian
# is scaled across the array. For a 64-pixel array, the 1/e points occur at approximately
# ±width×√2/2×63×64 ≈ ±width×0.7 pixels from center.

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

# Add reference lines and ticks
hlines!(ax, [1 / √ℯ]; color=:red, linestyle=:dot, linewidth=2, alpha=0.8)
vlines!(ax, [-4, 4]; color=:blue, linestyle=:dash, alpha=0.7)
vlines!(ax, [-7.5, 7.5]; color=:orange, linestyle=:dash, alpha=0.7)
vlines!(ax, [-12.5, 12.5]; color=:green, linestyle=:dash, alpha=0.7)

ax.yticks = ([1 / √ℯ], ["1/√e"])
ax.xticks = ([-12.5, -7.5, -4, 4, 7.5, 12.5], ["-12.5", "-7.5", "-4", "4", "7.5", "12.5"])

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

# ## Practical Application: Windowed Fourier Analysis
#
# Gaussian windows are commonly used to localize Fourier transforms. Here's an example
# showing how different window sizes affect spectral analysis:

using FFTW

## Increase resolution for better spectral analysis
N = 256
x = 1:N
# Add a slowly varying quadratic background to illustrate edge effects
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

# Compute FFT spectra
fft_original = abs.(fft(test_signal))
fft_narrow = abs.(fft(windowed_narrow))
fft_medium = abs.(fft(windowed_medium))
fft_wide = abs.(fft(windowed_wide))

fft_original ./= sum(abs2.(fft_original))
fft_narrow ./= sum(abs2.(fft_narrow))
fft_medium ./= sum(abs2.(fft_medium))
fft_wide ./= sum(abs2.(fft_wide))

freq = fftfreq(N, 1.0)

fig = Figure(; size=(1000, 600))

## Time domain
ax1 = Axis(fig[1, 1]; title="Windowed Signals", xlabel="Sample", ylabel="Amplitude")
lines!(ax1, x, test_signal; label="Original", linewidth=2)
lines!(ax1, x, windowed_narrow; label="Narrow window", linewidth=2)
lines!(ax1, x, windowed_medium; label="Medium window", linewidth=2)
lines!(ax1, x, windowed_wide; label="Wide window", linewidth=2)
axislegend(ax1; position=:rt)

## Frequency domain

ax2 = Axis(fig[2, 1]; title="FFT Spectra", xlabel="Frequency (in 1/N)", ylabel="Magnitude")
lines!(ax2, freq[1:div(N, 2)], fft_original[1:div(N, 2)]; label="Original", linewidth=2)
lines!(ax2, freq[1:div(N, 2)], fft_narrow[1:div(N, 2)]; label="Narrow window", linewidth=2)
lines!(ax2, freq[1:div(N, 2)], fft_medium[1:div(N, 2)]; label="Medium window", linewidth=2)
lines!(ax2, freq[1:div(N, 2)], fft_wide[1:div(N, 2)]; label="Wide window", linewidth=2)
f1_freq = f1 / N
f2_freq = f2 / N
## Add vertical lines for ground truth frequencies
vlines!(ax2, [f1_freq, f2_freq]; linewidth=2, alpha=0.7)
ax2.xticks = ([f1_freq, f2_freq], ["7.3", "3.7"])
axislegend(ax2; position=:rt)

fig


# ## Summary
#
# The `GaussianWindow` type in PhaseUtils provides flexible windowing capabilities:
#
# - **Symmetric windows**: `GaussianWindow(width)` for equal width in all dimensions
# - **Asymmetric windows**: `GaussianWindow((width_x, width_y, ...))` for dimension-specific widths
# - **Custom centering**: `center` parameter for off-center windows
# - **Backward compatibility**: Existing code continues to work unchanged
#
# These windows are particularly useful for:
# - Localizing Fourier analysis
# - Reducing edge effects in image processing
# - Creating smooth transitions in phase analysis
# - Controlling spectral leakage in signal processing
#
# The implementation maintains high performance while providing intuitive, flexible control
# over window shape and positioning.
