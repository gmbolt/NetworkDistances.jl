using Distances
export pre_compute, Precomputed

"""
    Precomputed{T<:SemiMetric,S}

Distance metric which looks up a cache of stored distances.

# Fields
- `d::T`: The underlying semi-metric for which distances were pre-computed.
- `cache::Dict{Tuple{Union{Nothing,S},Union{Nothing,S}},Float64}`: A dictionary storing the pre-computed distances.
"""
struct Precomputed{T<:SemiMetric,S} <: SemiMetric
    d::T
    cache::Dict{Tuple{Union{Nothing,S},Union{Nothing,S}},Float64}
end

function Base.show(io::IO, d_pre::Precomputed{T,S}) where {T<:SemiMetric,S}
    b = IOBuffer()
    show(b, d_pre.d)
    d_str = String(take!(b))
    print(io, "Precomputed{$(d_str),$(S)}")
end

"""
    (d::Precomputed{T,S})(args...) where {T<:SemiMetric,S}

Looks up and returns a pre-computed distance from the cache.

# Arguments
- `args...`: The arguments (typically two elements) for which to retrieve the distance.

# Returns
- `Float64`: The pre-computed distance.
"""
function (d::Precomputed{T,S})(
    args...
) where {T<:SemiMetric,S}
    return d.cache[args]
end

"""
    Precomputed(d::SemiMetric, data::Vector{S}) where {S}

Constructor for `Precomputed` that takes a semi-metric and a collection of values to pre-compute distances for.

# Arguments
- `d::SemiMetric`: The semi-metric to use for pre-computation.
- `data::Vector{S}`: A vector of values for which to pre-compute distances.

# Returns
- `Precomputed{typeof(d),S}`: A `Precomputed` distance object with the cache populated.
"""
function Precomputed(d::SemiMetric, data::Vector{S}) where {S}
    TypeVal = Union{S,Nothing} # We include distances to nothing.
    cache = Dict{Tuple{TypeVal,TypeVal},Float64}()
    for (a, b) in Base.Iterators.product([data; nothing], [data; nothing])
        cache[(a, b)] = d(a, b)
    end
    return Precomputed(d, cache)
end

"""
    pre_compute(d::SemiMetric, data::Vector{T}) where {T}

Convenience constructor method which takes a distance and a collection of values, and returns a `Precomputed` distance object.

# Arguments
- `d::SemiMetric`: The semi-metric to use for pre-computation.
- `data::Vector{T}`: A vector of values for which to pre-compute distances.

# Returns
- `Precomputed{typeof(d),T}`: A `Precomputed` distance object with the cache populated.
"""
function pre_compute(d::SemiMetric, data::Vector{T}) where {T}

    TypeVal = Union{T,Nothing}
    cache = Dict{Tuple{TypeVal,TypeVal},Float64}()
    for (a, b) in Base.Iterators.product([data; nothing], [data; nothing])
        cache[(a, b)] = d(a, b)
    end
    return Precomputed(d, cache)

end


