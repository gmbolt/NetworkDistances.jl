abstract type MatchingOptimiser end

"""
    ContinousRelaxation

Struct representing the continuous relaxation matching optimiser.
"""
struct ContinousRelaxation <: MatchingOptimiser end

"""
    HungarianAlgorithm

Struct representing the Hungarian algorithm matching optimiser.
"""
struct HungarianAlgorithm <: MatchingOptimiser end

"""
    eval_distance(optimiser::ContinousRelaxation, cost_matrix::AbstractMatrix)::Float64

Evaluates the distance using continuous relaxation.

# Arguments
- `optimiser::ContinousRelaxation`: The continuous relaxation optimiser.
- `cost_matrix::AbstractMatrix`: The cost matrix.

# Returns
- `Float64`: The evaluated distance.
"""
function eval_distance(optimiser::ContinousRelaxation, cost_matrix::AbstractMatrix)::Float64
    x = ones(size(cost_matrix, 1))
    return emd2(
        x, x, cost_matrix
    )
end

"""
    eval_distance(optimiser::HungarianAlgorithm, cost_matrix::AbstractMatrix)::Float64

Evaluates the distance using the Hungarian algorithm.

# Arguments
- `optimiser::HungarianAlgorithm`: The Hungarian algorithm optimiser.
- `cost_matrix::AbstractMatrix`: The cost matrix.

# Returns
- `Float64`: The evaluated distance.
"""
function eval_distance(optimiser::HungarianAlgorithm, cost_matrix::AbstractMatrix)::Float64
    _, cost = hungarian(cost_matrix)
    return cost
end

"""
    get_optimiser_instance(s::String)

Returns an instance of a `MatchingOptimiser` based on the provided string.

# Arguments
- `s::String`: The name of the optimiser ("hungarian" or "continuous_relaxation").

# Returns
- `MatchingOptimiser`: An instance of the specified optimiser.

# Throws
- `NotImplementedError`: If an unrecognised optimiser specification is provided.
"""
function get_optimiser_instance(s::String)
    optimiser = if (s == "hungarian")
        HungarianAlgorithm()
    elseif (s == "continuous_relaxation")
        ContinousRelaxation()
    else
        throw(NotImplementedError("Un-recognised optimiser specification. "))
    end
    return optimiser
end


