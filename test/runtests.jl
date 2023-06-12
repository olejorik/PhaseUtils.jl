using PhaseUtils
using Test

@testset "PhaseUtils.jl" begin
    # Write your tests here.

    aaa = reshape(1:9, (3, 3))

    ## Phase differentiations
    @test all(
        PhaseUtils._calculate_gradient(aaa, PhaseUtils.FiniteDifferencesCyclic()) .==  
        [[-6 3 3; -6 3 3; -6 3 3], [-2 -2 -2; 1 1 1; 1 1 1] ])
end

@testset "PhaseUnwrapping" begin
    s1,s2,m = 100,140,25
    wedge = [0.5i for i in 1:s1, j in 1:s2]
    wedge_ap  = [ m<i<s1-m && 25<j<s2-m for i in 1:s1, j in 1:s2]
    wedge .*= wedge_ap
    sol = unwrap_LS(phwrap(wedge), wedge_ap, restore_piston = true)
    @test all(phwrap(sol .* wedge_ap) .≈ phwrap(wedge .* wedge_ap))
end
