using Test
using PhaseUtils

@testset "axes_and_tilts: types and helpers" begin
    # Construct tilts
    t1 = TiltCentered([0.1, 0.2, -0.3])
    @test isa(t1, Tilt)
    @test sigma(t1) ≈ 0.1
    @test tau(t1) == [0.2, -0.3]
    @test tau(t1, 1) ≈ 0.2

    setsigma!(t1, 1.5)
    @test sigma(t1) ≈ 1.5
    settau!(t1, [-2.0, 3.0])
    @test all(tau(t1) .≈ [-2.0, 3.0])

    # FreeTilt behaves the same
    t2 = FreeTilt([0.0, 1.0, 2.0])
    @test sigma(t2) ≈ 0.0
    @test tau(t2, 2) ≈ 2.0

    # apply over a point
    @test apply(t1, [0.5, -0.5]) ≈ sigma(t1) + sum(tau(t1) .* [0.5, -0.5])

    # materialize over dims in Fourier-freq coordinates
    A = materialize(t1, (8, 10))
    @test size(A) == (8, 10)

    # axes policies
    dims = (6, 5)
    axF = FourierAxes()(dims)
    axD = DataAxes()(dims)
    axC = DataAxesCentered()(dims)

    @test length(axF) == 2 && length(axD) == 2 && length(axC) == 2
    @test length(axF[1]) == 6 && length(axF[2]) == 5
    @test axD[1] == 1:6 && axD[2] == 1:5

    # materialize with provided axes vector
    A2 = materialize(t1, axF)
    @test size(A2) == dims

    # sanity: materialize with axes equals materialize with dims for FourierAxes
    τ = tau(t1)
    A3 = [sigma(t1) + τ[1] * x + τ[2] * y for (x, y) in Iterators.product(axF[1], axF[2])]
    @test A2 ≈ A3

    # ProvidedAxes: wrap explicit axes and verify
    pax = ProvidedAxes(axF[1], axF[2])
    # Use the dims that correspond to the provided axes (NOT size(A) which is (8,10))
    A4 = materialize(t1, pax(dims))
    @test A4 ≈ A2
    # mismatch dims should throw
    @test_throws ArgumentError pax((length(axF[1]) + 1, length(axF[2])))
end
