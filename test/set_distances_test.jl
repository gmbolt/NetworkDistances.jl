@testset "Set distances tests" begin
    if isdefined(NetworkDistances, :set_distance)
        S1 = [1, 2, 3]
        S2 = [1, 2, 3]
        S3 = [4, 5, 6]
        @test set_distance(S1, S2) == 0.0
        @test set_distance(S1, S3) > 0.0
    else
        @info "set_distance not defined in NetworkDistances; skipping set distances tests"
    end
end 