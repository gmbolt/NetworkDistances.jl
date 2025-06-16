using Distances, StatsBase, ProgressMeter
export pairwise_inbounds, pairwise_inbounds!, progress_pairwise

"""
    Distances.pairwise!(A::AbstractArray, metric::Metric, a::Vector{T}, b::Vector{T})

Computes the distance matrix between elements of two vectors `a` and `b` using the given `metric`,
and stores the result in-place in `A`.

# Arguments
- `A::AbstractArray`: The output array to store the distance matrix.
- `metric::Metric`: The distance metric to use.
- `a::Vector{T}`: The first vector of elements.
- `b::Vector{T}`: The second vector of elements.
"""
function Distances.pairwise!(
    A::AbstractArray,
    metric::Metric,
    a::Vector{T} where {T},
    b::Vector{T} where {T}
)

    for j in eachindex(b)
        for i in eachindex(a)
            A[i, j] = metric(a[i], b[j])
        end
    end
end

"""
    pairwise_inbounds!(A::AbstractArray, metric::Metric, a::Vector{T}, b::Vector{T}) where {T}

Computes the distance matrix between elements of two vectors `a` and `b` using the given `metric`,
and stores the result in-place in `A`, with `@inbounds` optimization.

# Arguments
- `A::AbstractArray`: The output array to store the distance matrix.
- `metric::Metric`: The distance metric to use.
- `a::Vector{T}`: The first vector of elements.
- `b::Vector{T}`: The second vector of elements.
"""
function pairwise_inbounds!(
    A::AbstractArray,
    metric::Metric,
    a::Vector{T},
    b::Vector{T}
) where {T}
    @inbounds begin
        for j in eachindex(b)
            for i in eachindex(a)
                A[i, j] = metric(a[i], b[j])
            end
        end
    end
end

"""
    Distances.pairwise!(A::SubArray, metric::Metric, a::Vector{T}, b::Vector{T}) where {T}

Computes the distance matrix between elements of two vectors `a` and `b` using the given `metric`,
and stores the result in-place in `A` (for SubArray).

# Arguments
- `A::SubArray`: The output SubArray to store the distance matrix.
- `metric::Metric`: The distance metric to use.
- `a::Vector{T}`: The first vector of elements.
- `b::Vector{T}`: The second vector of elements.
"""
function Distances.pairwise!(
    A::SubArray,
    metric::Metric,
    a::Vector{T} where {T},
    b::Vector{T} where {T}
)

    for j in eachindex(b)
        for i in eachindex(a)
            A[i, j] = metric(a[i], b[j])
        end
    end
end


"""
    Distances.pairwise!(A::AbstractMatrix, metric::Metric, a::Vector{T}, b::Vector{T}) where {T}

Computes the distance matrix between elements of two vectors `a` and `b` using the given `metric`,
and stores the result in-place in `A` (for AbstractMatrix).

# Arguments
- `A::AbstractMatrix`: The output AbstractMatrix to store the distance matrix.
- `metric::Metric`: The distance metric to use.
- `a::Vector{T}`: The first vector of elements.
- `b::Vector{T}`: The second vector of elements.
"""
function Distances.pairwise!(
    A::AbstractMatrix,
    metric::Metric,
    a::Vector{T} where {T},
    b::Vector{T} where {T}
)

    for j in eachindex(b)
        for i in eachindex(a)
            A[i, j] = metric(a[i], b[j])
        end
    end
end


"""
    Distances.pairwise(metric::SemiMetric, a::Vector{T}, b::Vector{T}) where {T}

Computes the distance matrix between elements of two vectors `a` and `b` using the given `metric`.
This is a custom extension of the function in the Distances.jl package to allow vectors of general type.
The function in Distances.jl is designed for univariate/multivariate data and so takes
as input either vectors or matrices (data points as rows).

# Arguments
- `metric::SemiMetric`: The distance metric to use.
- `a::Vector{T}`: The first vector of elements.
- `b::Vector{T}`: The second vector of elements.

# Returns
- `Matrix{Float64}`: The computed distance matrix.
"""
function Distances.pairwise(
    metric::SemiMetric,
    a::Vector{T},
    b::Vector{T}
) where {T}
    D = Array{Float64,2}(undef, length(a), length(b))
    for j in eachindex(b)
        for i in eachindex(a)
            D[i, j] = metric(a[i], b[j])
        end
    end
    return D
end

"""
    pairwise_inbounds(metric::SemiMetric, a::Vector{T}, b::Vector{T}) where {T}

Computes the distance matrix between elements of two vectors `a` and `b` using the given `metric`,
with `@inbounds` optimization.

# Arguments
- `metric::SemiMetric`: The distance metric to use.
- `a::Vector{T}`: The first vector of elements.
- `b::Vector{T}`: The second vector of elements.

# Returns
- `Matrix{Float64}`: The computed distance matrix.
"""
function pairwise_inbounds(
    metric::SemiMetric,
    a::Vector{T},
    b::Vector{T}
) where {T}
    D = Array{Float64,2}(undef, length(a), length(b))
    @inbounds begin
        for j in eachindex(b)
            for i in eachindex(a)
                D[i, j] = metric(a[i], b[j])
            end
        end
    end
    return D
end

"""
    progress_pairwise(d::SemiMetric, a::Vector{T}) where {T}

Computes the pairwise distance matrix for a single vector `a` using the given `SemiMetric` `d`,
and displays a progress bar.

# Arguments
- `d::SemiMetric`: The distance metric to use.
- `a::Vector{T}`: The vector of elements.

# Returns
- `Matrix{Float64}`: The symmetric distance matrix.
"""
function progress_pairwise(
    d::SemiMetric,
    a::Vector{T}
) where {T}

    D = zeros(length(a), length(a))
    iter = Progress(Int(length(a) * (length(a) - 1) / 2), 1)
    for j in eachindex(a)
        for i in 1:(j-1)
            D[i, j] = d(a[i], a[j])
            next!(iter)
        end
    end
    D += D'
    return D

end


