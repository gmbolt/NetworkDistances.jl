using Test, NetworkDistances


@testset "NetworkDistances.jl" begin

    println(readdir())
    println(readdir("distances/"))
    include("matching_distances_test.jl")
    include("distances/multisets_new/matching_distances_test.jl")



end
