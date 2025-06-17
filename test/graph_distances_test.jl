using Test
using NetworkDistances

@testset "Graph Distances" begin
    # Define two simple adjacency matrices
    A = [0 1; 1 0]
    B = [0 1; 0 0]

    # Check self-distance is zero
    @test hamming_dist(A, A) == 0.0
    @test jaccard_dist(A, A) == 0.0
    @test diffusion_dist(A, A, 0.5) == 0.0

    # Check symmetry for non-identical matrices
    d1 = hamming_dist(A, B)
    d2 = hamming_dist(B, A)
    @test d1 == d2

    d1 = jaccard_dist(A, B)
    d2 = jaccard_dist(B, A)
    @test d1 == d2

    d1 = diffusion_dist(A, B, 0.5)
    d2 = diffusion_dist(B, A, 0.5)
    @test isapprox(d1, d2, atol=1e-10)
end