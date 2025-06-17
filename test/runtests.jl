using Test, NetworkDistances

@testset "NetworkDistances.jl" begin
    include("matching_distances_test.jl")
    include("graph_distances_test.jl")
    include("path_distances_test.jl")
    include("multiset_emd_test.jl")
    include("sequence_distances_test.jl")
    include("set_distances_test.jl")
    include("utils_test.jl")
end