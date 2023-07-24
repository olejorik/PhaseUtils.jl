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
    end

    @testset "PhaseUnwrapping" begin
        s1, s2, m = 140, 100, 25
        wedge = [0.5i for i in 1:s1, j in 1:s2]
        wedge_ap = [m < i < s1 - m && 25 < j < s2 - m for i in 1:s1, j in 1:s2]
        wedge .*= wedge_ap
        sol = unwrap_LS(phwrap(wedge), wedge_ap; restore_piston=true)
        @test all(phwrap(sol .* wedge_ap) .≈ phwrap(wedge .* wedge_ap))

        # Example with a circular aperture
        ap = zeros(s1, s2)
        y = range(-1.1 * s1 / s2, 1.1 * s1 / s2, s1)
        x = range(-1.1, 1.1, s2)
        ap[[x .^ 2 + y .^ 2 .<= 1 for y in y, x in x]] .= 1
        phaseGT = [-x^3 + 3x .^ 2 + y .^ 2 - 10y for y in y, x in x]
        phase = phwrap(phaseGT) .* ap
        sol = unwrap_LS(phase, ap; restore_piston=true)
        @test all(phwrap(phaseGT .- sol) .* ap .+ 1 .≈ 1.0)  # the difference can be 2pi*k, so we wrap it
        # relatively to 1., the error should be negligible.
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
end
