using StatsBase, Distances, Hungarian, Printf, LogExpFunctions

export MatchingDistance, MatchDist
export FastMatchingDistance, FastMatchDist
export FixedPenaltyMatchingDistance, FixPenMatchDist
export MinDistMatchingDistance, MinDistMatchDist
export AvgSizeMatchingDistance, AvgSizeMatchDist
export get_cost_matrix_dynamic, get_cost_matrix_fixed

export FixedPenalty, SizePenalty, DistancePenalty, ParametricPenalty

abstract type AbstractMatchingDistance <: SemiMetric end
struct NotImplementedError <: Exception
    msg::String
end

"""
    get_cost_matrix_dynamic(d::AbstractMatchingDistance, X, Y)

Return cost matrix for evaluating matching-based distance between `X` and `Y`. 

This is done "dynamically" in the sense that the first axis of the output matrix 
will always correspond to the shorter of the two passed vectors `X` and `Y`.
"""
function get_cost_matrix_dynamic(
    d::AbstractMatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}
    return NotImplementedError("Method get_cost_matrix_dynamic() not implemented for this distance.")
end

"""
    get_cost_matrix_fixed(d::AbstractMatchingDistance, X, Y)

Return cost matrix for evaluating matching-based distance between `X` and `Y`. 

This is "fixed" in the sense the first axis of output matrix always corresponds to `X` and the second to `Y`.
"""
function get_cost_matrix_fixed(
    d::AbstractMatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}
    return NotImplementedError("Method get_cost_matrix_fixed() not implemented for this distance.")
end

abstract type MatchingOptimiser end
struct ContinousRelaxation <: MatchingOptimiser end
struct HungarianAlgorithm <: MatchingOptimiser end

function eval_distance(optimiser::ContinousRelaxation, cost_matrix::AbstractMatrix)::Float64
    x = ones(size(cost_matrix, 1))
    return PythonOT.emd2(
        x, x, cost_matrix
    )
end

function eval_distance(optimiser::HungarianAlgorithm, cost_matrix::AbstractMatrix)::Float64
    assignment, cost = hungarian(cost_matrix)
    return cost
end

abstract type CompleteMatchingDistance <: AbstractMatchingDistance end

# We will default to continuous relaxation for evaluation of distances
get_optimiser(d::CompleteMatchingDistance) = ContinousRelaxation()

abstract type PenaltyFunction end

struct FixedPenalty <: PenaltyFunction
    rho::Float64
end

(penalty::FixedPenalty)(x::Vector{T}) where {T} = penalty.rho

Base.show(io::IO, p::S) where {S<:FixedPenalty} = print(io, "$(S)(ρ=$(p.rho))")

struct SizePenalty <: PenaltyFunction end

(penalty::SizePenalty)(x::Vector{T}) where {T} = length(x)

Base.show(io::IO, p::S) where {S<:SizePenalty} = print(io, "$(S)")

struct DistancePenalty{T<:SemiMetric} <: PenaltyFunction
    d::T
end

(penalty::DistancePenalty)(x::Vector{T}) where {T} = penalty.d(x, nothing)

Base.show(io::IO, p::DistancePenalty{S}) where {S<:SemiMetric} = print(
    io, "DistancePenalty{$(S)}"
)

struct ParametricPenalty <: PenaltyFunction
    loc::Float64
    scale::Float64
    interc::Float64
    function ParametricPenalty(loc::Real, scale::Real; interc::Real=0.0)
        new(loc, scale, interc)
    end
end

function (penalty::ParametricPenalty)(x::Vector{T}) where {T}
    x_len = length(x)
    loc, scale, interc = (penalty.loc, penalty.scale, penalty.interc)
    return (
        scale * x_len
        +
        scale * xlogx(loc)
        -
        scale * loc * (log(x_len) + 1)
        +
        interc
    )
end

Base.show(io::IO, p::ParametricPenalty) = print(
    io, "ParametricPenalty(loc=$(p.loc),scale=$(p.scale))"
)


struct MatchingDistance{T<:SemiMetric,S<:PenaltyFunction,R<:MatchingOptimiser} <: CompleteMatchingDistance
    ground_dist::T
    penalty::S
    optimiser::R
end

# Optimier query function
get_optimiser(d::MatchingDistance) = d.optimiser
# Constructors
MatchingDistance(d::SemiMetric, penalty::PenaltyFunction; optimiser::MatchingOptimiser=ContinousRelaxation()) = MatchingDistance(d, penalty, optimiser)
MatchingDistance(d::SemiMetric; optimiser::MatchingOptimiser=ContinousRelaxation()) = MatchingDistance(d, DistancePenalty(d), optimiser)

const MatchDist = MatchingDistance

function get_cost_matrix_dynamic(
    d::MatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}

    if length(X) < length(Y)
        return get_cost_matrix_dynamic(d, Y, X)
    elseif length(X) == length(Y)
        return pairwise_inbounds(d.ground_dist, X, Y)
    else
        C = pairwise_inbounds(d.ground_dist, X, Y)
        null_dists = [d.penalty(p) for p in X]
        size_diff = length(X) - length(Y)
        return [C [x for x ∈ null_dists, j = 1:size_diff]]
    end
end

function get_cost_matrix_fixed(
    d::MatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}

    if length(X) < length(Y)
        C = pairwise_inbounds(d.ground_dist, X, Y)
        pentalty_vec = [d.ground_dist(nothing, p) for p in Y]
        size_diff = length(Y) - length(X)
        return [C; [y for i = 1:size_diff, y ∈ pentalty_vec]]
    elseif length(X) == length(Y)
        return pairwise_inbounds(d.ground_dist, X, Y)
    else
        C = pairwise_inbounds(d.ground_dist, X, Y)
        penalty_vec = [d.penalty(p) for p in X]
        size_diff = length(X) - length(Y)
        return [C [x for x ∈ penalty_vec, j = 1:size_diff]]
    end

end


struct FastMatchingDistance{T<:SemiMetric,S<:PenaltyFunction,R<:MatchingOptimiser} <: CompleteMatchingDistance
    ground_dist::T
    penalty::S
    optimiser::R
    C::Matrix{Float64}
    function FastMatchingDistance(ground_dist::T, penalty::S, optimiser::R, K::Int) where {T<:SemiMetric,S<:PenaltyFunction,R<:MatchingOptimiser}
        new{T,S,R}(ground_dist, penalty, optimiser, zeros(K, K))
    end
end


const FastMatchDist = FastMatchingDistance
FastMatchingDistance(d::SemiMetric, penalty::PenaltyFunction, K::Int; optimiser::MatchingOptimiser=ContinousRelaxation()) = FastMatchingDistance(d, penalty, optimiser, K)
FastMatchingDistance(d::SemiMetric, K::Int; optimiser::MatchingOptimiser=ContinousRelaxation()) = FastMatchingDistance(d, DistancePenalty(d), optimiser, K)

Base.show(io::IO, d::FastMatchingDistance{T,S,R}) where {T<:SemiMetric,S<:PenaltyFunction,R<:MatchingOptimiser} = print(io, "FastMatchingDistance{$(T),$(S),$(R),$(size(d.C,1))}")

function get_cost_matrix_dynamic(
    d::FastMatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}

    N, M = (length(X), length(Y))
    if N < M
        return get_cost_matrix_dynamic(d, Y, X)
    elseif N == M
        C = view(d.C, 1:N, 1:N)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
    else
        C = view(d.C, 1:N, 1:N)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
        for i in 1:N
            C[i, M+1] = d.penalty(X[i])
        end
        for j in (M+2):N
            for i in 1:N
                C[i, j] = C[i, j-1]
            end
        end
    end
    return C
end

function get_cost_matrix_fixed(
    d::FastMatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}

    N, M = (length(X), length(Y))
    if N < M
        C = view(d.C, 1:M, 1:M)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
        for j in 1:M
            C[N+1, j] = d.penalty(Y[j])
        end
        for i in (N+2):N
            for j in 1:M
                C[i, j] = C[i-1, j]
            end
        end
    elseif N == M
        C = view(d.C, 1:N, 1:N)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
    else
        C = view(d.C, 1:N, 1:N)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
        for i in 1:N
            C[i, M+1] = d.penalty(X[i])
        end
        for j in (M+2):M
            for i in 1:N
                C[i, j] = C[i, j-1]
            end
        end
    end
    return C

end

function (d::Union{MatchDist,FastMatchDist})(X::Nothing, Y::Vector{T})::Float64 where {T}
    return sum(p -> d.penalty(p), Y)
end
function (d::Union{MatchDist,FastMatchDist})(X::Vector{T}, Y::Nothing)::Float64 where {T}
    return d(Y, X)
end
function (d::Union{MatchDist,FastMatchDist})(X::Nothing, Y::Nothing)::Float64
    return 0.0
end

# Legacy values for backwards compatability
const FixedPenaltyMatchingDistance{T} = MatchingDistance{T,FixedPenalty} where {T<:SemiMetric}
const FixPenMatchDist{T} = FixedPenaltyMatchingDistance{T} where {T<:SemiMetric}
FixedPenaltyMatchingDistance(d::SemiMetric, penalty::Real) = MatchDist(d, FixedPenalty(penalty))
FixPenMatchDist(args...) = FixedPenaltyMatchingDistance(args...)

struct AvgSizeMatchingDistance{T<:SemiMetric} <: CompleteMatchingDistance
    ground_dist::T
    penalty::Float64
end

const AvgSizeMatchDist = AvgSizeMatchingDistance

Base.show(io::IO, d::AvgSizeMatchDist) = print(io, "$(typeof(d))(ρ=$(d.penalty))")

function get_cost_matrix_dynamic(
    d::AvgSizeMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    if N < M
        return get_cost_matrix_dynamic(d, Y, X)
    elseif N == M
        return pairwise_inbounds(d_g, X, Y)
    else
        C = pairwise_inbounds(d_g, X, Y)
        mean_size = mean(y -> d_g(y, nothing), Y) # Find mean size in smaller set
        pentalty_vec = [abs(d_g(x, nothing) - mean_size) + d.penalty for x in X]
        size_diff = (N - M)
        return [C [x for x ∈ pentalty_vec, j = 1:size_diff]]
    end
end

function get_cost_matrix_fixed(
    d::AvgSizeMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    C = pairwise_inbounds(d_g, X, Y)
    if N < M
        mean_size = mean(x -> d_g(x, nothing), X)
        penalty_vec = [abs(d_g(y, nothing) - mean_size) + d.penalty for y in Y]
        size_diff = (M - N)
        return [C; [y for i in 1:size_diff, y ∈ penalty_vec]]
    elseif N > M
        mean_size = mean(y -> d_g(y, nothing), Y)
        penalty_vec = [abs(d_g(x, nothing) - mean_size) + d.penalty for x in X]
        size_diff = (N - M)
        return [C [x for x ∈ penalty_vec, i in 1:size_diff]]
    else
        return pairwise_inbounds(d_g, X, Y)
    end

end

function (d::AvgSizeMatchDist)(X::Nothing, Y::Vector{T})::Float64 where {T}
    return (d.penalty * length(Y)) + sum(x -> dist.ground_dist(x, nothing), Y)
end
function (d::AvgSizeMatchDist)(X::Vector{T}, Y::Nothing)::Float64 where {T}
    d(Y, X)
end
function (d::AvgSizeMatchDist)(X::Nothing, Y::Nothing)::Float64
    return 0.0
end

struct MinDistMatchingDistance{T<:SemiMetric} <: CompleteMatchingDistance
    ground_dist::T
    penalty::Float64
end

const MinDistMatchDist{T} = MinDistMatchingDistance where {T}

Base.show(io::IO, d::MinDistMatchDist) = print(io, "$(typeof(d))(ρ=$(d.penalty))")

function get_cost_matrix_dynamic(
    d::MinDistMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    if N < M
        return get_cost_matrix_dynamic(d, Y, X)
    elseif N == M
        return pairwise_inbounds(d_g, X, Y)
    else
        C = pairwise_inbounds(d_g, X, Y)
        penalty_vec = map(minimum, eachrow(C))
        size_diff = (N - M)
        return [C [x for x ∈ penalty_vec, j = 1:size_diff]]
    end
end

function get_cost_matrix_fixed(
    d::MinDistMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    C = pairwise_inbounds(d_g, X, Y)
    if N < M
        penalty_vec = map(minimum, eachcol(C))
        size_diff = (M - N)
        return [C; [y for i in 1:size_diff, y ∈ penalty_vec]]
    elseif N > M
        penalty_vec = map(minimum, eachrow(C))
        size_diff = (N - M)
        return [C [x for x ∈ penalty_vec, i in 1:size_diff]]
    else
        return pairwise_inbounds(d_g, X, Y)
    end

end
