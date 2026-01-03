@testset "matching_distances" begin
    V = 1:10
    n = 30
    x = [rand(V, rand(4:7)) for i in 1:n]
    y = [rand(V, rand(4:7)) for i in 1:2n]
    d1 = MatchingDistance(LCS())
    d2 = FastMatchDist(LCS(), 2n)
    d3 = FixPenMatchDist(LCS(), 1.0)
    # Test that d1 and d2 produce similar results for the test data
    val1 = d1(x, y)
    val2 = d2(x, y)
    @test isapprox(val1, val2, atol=1e-8)
    # Check that distances are zero for identical inputs
    @test d1(x, x) == 0.0
    @test d2(x, x) == 0.0
    @test d3(nothing, nothing) == 0.0
    # Check symmetry
    @test isapprox(d1(x, y), d1(y, x), atol=1e-8)
    @test isapprox(d2(x, y), d2(y, x), atol=1e-8)
    @test isapprox(d3(x, y), d3(y, x), atol=1e-8)
end
