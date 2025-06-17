using Test
using NetworkDistances

@testset "Multiset EMD tests" begin
    # Use AbsoluteDiff as ground distance for EMD.
    ground = AbsoluteDiff()
    d = EMD(ground)
    X = [1, 1, 2, 2]
    Y = [1, 2, 2, 3]

    dist = d(X, Y)
    @test dist >= 0.0

    # Test get_info for EMD: should return (keysA, keysB, transport plan)
    a_keys, b_keys, transPlan = get_info(d, X, Y)
    @test length(a_keys) > 0
    @test length(b_keys) > 0
    @test size(transPlan, 1) == length(a_keys)
    @test size(transPlan, 2) == length(b_keys)
end