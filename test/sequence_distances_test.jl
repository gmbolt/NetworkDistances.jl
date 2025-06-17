using Test
using NetworkDistances

@testset "Sequences distances tests" begin
    # DTW Distance tests
    series1 = [1.0, 2.0, 3.0]
    series2 = [1.0, 2.0, 3.0]
    series3 = [2.0, 3.0, 4.0]

    @test dtw_distance(series1, series2) == 0.0
    d1 = dtw_distance(series1, series3)
    @test d1 > 0.0
    @test dtw_distance(series1, series3) == dtw_distance(series3, series1)

    # Edit Distance tests (for strings)
    s1 = "kitten"
    s2 = "kitten"
    s3 = "sitting"
    @test edit_distance(s1, s2) == 0
    @test edit_distance(s1, s3) == 3

    # Edit Distance tests (for arrays)
    arr1 = [1, 2, 3, 4]
    arr2 = [1, 3, 4]
    d2 = edit_distance(arr1, arr2)
    @test d2 > 0
end 