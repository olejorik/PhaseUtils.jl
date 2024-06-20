using PhaseUtils
using Test

@testset "PhaseUtils.jl" begin
    # Write your tests here.

    @testset "utils" begin
        ttt = ap2mask([0.5, 0])
        @test all([ttt[1] == 1, isnan(ttt[2])])
        @test all(mask2ap(ttt) .== [1, 0])
    end

    @testset "differentiations" begin
        aaa = reshape(1:9, (3, 3))

        ## Phase differentiations
        @test all(
            PhaseUtils._calculate_gradient(aaa, PhaseUtils.FiniteDifferencesCyclic()) .==
            [[-6 3 3; -6 3 3; -6 3 3], [-2 -2 -2; 1 1 1; 1 1 1]],
        )

        ## Laplacian
        x = y = -1.2:0.2:1.2
        ap = @. x^2 + y'^2 <= 1
        f = @. (x^2 + y'^2 - 1) * ap
        l = PhaseUtils._calculate_Laplacian(f, PhaseUtils.FiniteDifferencesCyclic())
        @test all(l[1:2, 1:2] .== 0)
        @test all(l[6:8, 6:8] .≈ -0.16)
    end

    @testset "PhaseUnwrappingContour" begin
        # ϕ = 20 * cos.(range(0, 6π, 200)) .+ 5
        # w = phwrap(ϕ)
        # fig, ax, l = lines(ϕ; label="Ground truth", linewidth=5, color=(:gray, 0.5))
        # lines!(w; label="wrapped phase")
        # fi = itoh(w)
        # lines!(fi; label="Itoh algorithm")

        # fc = PhaseUtils.unwrap_contour_LS(w)
        # lines!(fc; label="Fourier algorithm")

        # axislegend(ax)
        # fig
        ϕ = 20 * cos.(range(0, 6π, 200))
        ϕ .-= (sum(ϕ) / length(ϕ))
        w = phwrap(ϕ)
        fc = PhaseUtils.unwrap_contour_LS(w; restore_piston=false)
        @test all(fc .≈ ϕ)
    end
    @testset "PhaseUnwrapping" begin
        s1, s2, m = 140, 100, 25
        wedge = [0.5i for i in 1:s1, j in 1:s2]
        wedge_ap = [m < i < s1 - m && 25 < j < s2 - m for i in 1:s1, j in 1:s2]
        wedge .*= wedge_ap
        wedge .-= sum(wedge) / sum(wedge_ap) #extract mean
        sol = unwrap_LS(phwrap(wedge), wedge_ap; restore_piston=true)
        sol .-= (sol[70, 50] - wedge[70, 50]) #extract 2kπ piston
        @test maskedrmse(sol, wedge, wedge_ap) / sum(wedge_ap) < 1e-16

        @test itoh(phwrap(1:100)) == 1:100

        # Example with a circular aperture
        ap = zeros(s1, s2)
        y = range(-1.1 * s1 / s2, 1.1 * s1 / s2, s1)
        x = range(-1.1, 1.1, s2)
        ap[[x .^ 2 + y .^ 2 .<= 1 for y in y, x in x]] .= 1
        phaseGT = [-x^3 + 3x .^ 2 + y .^ 2 - 10y for y in y, x in x]
        phase = phwrap(phaseGT) .* ap#extract mean
        sol = unwrap_LS(phase, ap; restore_piston=true)
        sol .-= sum(sol) / sum(ap)
        err = (sol .- phaseGT)
        err .-= sum(err .* ap) / sum(ap)
        @test maskedrmse(err, ap) / sum(ap) < 1e-16
    end

    @testset "Crop and Pad" begin
        a = reshape(1:15, 3, 5)
        @test PhaseUtils.croprange((2, 2), 2) == CartesianIndices((1:2, 1:2))
        @test PhaseUtils.croprange((2, 2), (1, 3)) == CartesianIndices((0:1, 2:3))
        @test crop(a, 1, (1, 3)) == [7;;]
        @test crop(a, (2, 2), (1, 3)) == [4 7;]
        @test crop(a, (2, 2), 2) == [1 4; 2 5]
        @test crop(a, 1) == [8;;]
        @test crop(a, 2) == [4 7; 5 8]
        @test crop(a, 3) == [4 7 10; 5 8 11; 6 9 12]
        @test crop(a, 2, 2) == [1 4; 2 5]
        @test crop(2, 2)(a) == [1 4; 2 5]
        @test crop(3)(a) == [4 7 10; 5 8 11; 6 9 12]
    end

    @testset "Vectorization" begin
        a = rand(1:10, (10, 10))
        circlemask!(a, 3)
        a[5, 10] = 1
        ind = findall(a .> 0)
        vals = a[ind]

        b = toArray(ind, vals; crop=false)
        @test all(a[1:8, :] .== b)

        b = toArray(ind, vals; crop=true)

        @test all(a[3:8, 3:end] .== b)

        b = toArray(ind, vals; shift=(0, 1))
        @test all(a[3:8, 3:end] .== b[:, 2:end])

        anan = Float64.(copy(a))
        anan[anan .== 0] .= NaN
        bnan = toArray(ind, vals; shift=(0, 1), defaultel=NaN)
        anansmall = anan[3:8, 2:end]
        @test all(anansmall[(!isnan).(anansmall)] .== bnan[(!isnan).(anansmall)]) &&
            all((isnan).(anansmall) .== (isnan).(bnan))
    end

    @testset "Thresholding" begin
        a = fill(1, (10, 10))
        circlemask!(a, 3)
        mask = copy(a)
        a .= a .* 10 .- 3

        alg = PhaseUtils.HardThreshold(5)
        b = PhaseUtils.threshold(a, alg)
        @test all(mask .* 7 .== b)

        b = PhaseUtils.threshold(a, PhaseUtils.HardThreshold(1.0))
        @test all(a .== b)


        b = PhaseUtils.threshold(a, PhaseUtils.SoftThreshold(5))
        @test all(mask .* 2 .== b)

        b = PhaseUtils.threshold(a, PhaseUtils.SoftThreshold(1.0))
        @test all(mask .* 8 .== b .+ 2)



    end
end
