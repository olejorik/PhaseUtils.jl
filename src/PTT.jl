# Simple utils for finding linear terms
#

function lineararray(size, kx, ky, k0=0)
    x = (1:size)'
    y = 1:size
    return kx * x .+ ky * y .+ k0
end
function lineararray(xrange::AbstractRange, yrange::AbstractRange, kx, ky, k0=0)
    return kx * xrange' .+ ky * yrange .+ k0
end

function lineararray(xrange::AbstractRange, yrange::AbstractRange, a::Vector, k0=0)
    return lineararray(xrange, yrange, a..., k0)
end

function findpiston(ϕ)
    return mean(filter(!isnan, ϕ))
end

function findpiston(ϕ, mask)
    return mean(phwrap.(filter(!isnan, ϕ .* mask)))
end

function findtiptilt(ϕ)
    sy, sx = size(ϕ)
    # dx = diff(ϕ; dims=2)
    # dy = diff(ϕ; dims=1)
    # kx = mean(phwrap.(filter(!isnan, dx)))
    # ky = mean(phwrap.(filter(!isnan, dy)))
    kx, ky = findtiptilt_coef(ϕ)
    tiptilt = lineararray(1:sx, 1:sy, kx, ky)
    return tiptilt
end

findtiptilt(ϕ, mask) = findtiptilt(ϕ .* mask)

function findtiptilt_coef(ϕ)
    dx = diff(ϕ; dims=2)
    dy = diff(ϕ; dims=1)
    kx = mean(phwrap.(filter(!isnan, dx)))
    ky = mean(phwrap.(filter(!isnan, dy)))
    return (kx, ky)
end

findtiptilt_coef(ϕ, mask) = findtiptilt_coef(ϕ .* mask)
