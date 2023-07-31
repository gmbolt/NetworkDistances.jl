using Test, NetworkDistances


@testset "NetworkDistances.jl" begin

    println(readdir())
    include("matching_distances_test.jl")
    include(joinpath(
        "distances",
        "multisets",
        "matching_distances_test.jl"
    ))



end
