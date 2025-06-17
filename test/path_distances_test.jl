using Test
using NetworkDistances

@testset "Paths distances tests" begin
    # Test sequences for LCS and FastLCS
    seq1 = [1, 2, 3, 4, 5]
    seq2 = [1, 3, 4, 5]
    seq3 = [5, 4, 3, 2, 1]
    tuple1 = ("a", "b", "c")
    tuple2 = ("a", "x", "c")

    # Test LCS (basic)
    lcs = LCS()
    @test lcs(seq1, seq1) == 0.0
    @test lcs(tuple1, tuple1) == 0.0
    @test lcs(seq1, seq2) == 1.0   # One element difference

    # Test FastLCS (optimized version) on vector input
    k = max(length(seq1), length(seq2))
    fastlcs = FastLCS(k)
    @test fastlcs(seq1, seq1) == 0.0
    @test fastlcs(seq1, seq2) == 1.0

    # Test behavior when one input is nothing
    @test lcs(nothing, seq1) == length(seq1)
    @test lcs(seq1, nothing) == length(seq1)
    @test lcs(nothing, nothing) == 0.0

    # Test get_info returns boolean arrays of proper length
    indx, indy = get_info(lcs, seq1, seq2)
    @test length(indx) == length(seq1)
    @test length(indy) == length(seq2)
    @test all(x -> (x == true || x == false), indx)
    @test all(x -> (x == true || x == false), indy)

    # Test LSP: identical sequences yield 0 distance; dissimilar ones yield > 0
    lsp = LSP()
    @test lsp(seq1, seq1) == 0.0
    @test lsp(seq1, seq3) > 0.0
end