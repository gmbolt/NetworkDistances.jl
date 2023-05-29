abstract type MatchingOptimiser end
struct ContinousRelaxation <: MatchingOptimiser end
struct HungarianAlgorithm <: MatchingOptimiser end

function eval_distance(optimiser::ContinousRelaxation, cost_matrix::AbstractMatrix)::Float64
    x = ones(size(cost_matrix, 1))
    return emd2(
        x, x, cost_matrix
    )
end

function eval_distance(optimiser::HungarianAlgorithm, cost_matrix::AbstractMatrix)::Float64
    _, cost = hungarian(cost_matrix)
    return cost
end

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