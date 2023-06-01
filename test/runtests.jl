using PhaseUtils
using Test

@testset "PhaseUtils.jl" begin
    # Write your tests here.

    aaa = reshape(1:9, (3, 3))

    ## Phase differentiations
    PhaseUtils._calculate_gradient(aaa, PhaseUtils.FiniteDifferencesCyclic()) ==
    ([-2 -2 -2; 1 1 1; 1 1 1], [-6 3 3; -6 3 3; -6 3 3])
end
