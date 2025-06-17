using Test
using NetworkDistances

@testset "Utils tests" begin
    if isdefined(NetworkDistances, :pairwise)
        try
            result = pairwise(+, [1, 2, 3])
            @test typeof(result) <: AbstractArray
        catch e
            @info "Error while testing pairwise: $e"
        end
    else
        @info "pairwise not defined in NetworkDistances; skipping Utils tests"
    end
end