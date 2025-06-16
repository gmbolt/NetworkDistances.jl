export FixedPenalty, SizePenalty, DistancePenalty, ParametricPenalty

abstract type PenaltyFunction end

"""
    FixedPenalty

Struct representing a fixed penalty function.

# Fields
- `rho::Float64`: The fixed penalty value.
"""
struct FixedPenalty <: PenaltyFunction
    rho::Float64
end

"""
    (penalty::FixedPenalty)(x::Vector{T}) where {T}

Applies the fixed penalty to any input.

# Arguments
- `x::Vector{T}`: The input vector (ignored).

# Returns
- `Float64`: The fixed penalty value `rho`.
"""
(penalty::FixedPenalty)(x::Vector{T}) where {T} = penalty.rho

Base.show(io::IO, p::S) where {S<:FixedPenalty} = print(io, "$(S)(ρ=$(p.rho))")

"""
    SizePenalty

Struct representing a penalty function based on the size (length) of the input.
"""
struct SizePenalty <: PenaltyFunction end

"""
    (penalty::SizePenalty)(x::Vector{T}) where {T}

Applies the size penalty, returning the length of the input vector.

# Arguments
- `x::Vector{T}`: The input vector.

# Returns
- `Int`: The length of the input vector.
"""
(penalty::SizePenalty)(x::Vector{T}) where {T} = length(x)

Base.show(io::IO, p::S) where {S<:SizePenalty} = print(io, "$(S)")

"""
    DistancePenalty{T<:SemiMetric}

Struct representing a penalty function based on a distance to `nothing`.

# Fields
- `d::T`: The semi-metric used to calculate the distance.
"""
struct DistancePenalty{T<:SemiMetric} <: PenaltyFunction
    d::T
end

"""
    (penalty::DistancePenalty)(x::Vector{T}) where {T}

Applies the distance penalty, returning the distance from `x` to `nothing` using the specified semi-metric.

# Arguments
- `x::Vector{T}`: The input vector.

# Returns
- `Float64`: The distance from `x` to `nothing`.
"""
(penalty::DistancePenalty)(x::Vector{T}) where {T} = penalty.d(x, nothing)

Base.show(io::IO, p::DistancePenalty{S}) where {S<:SemiMetric} = print(
    io, "DistancePenalty{$(S)}"
)

"""
    ParametricPenalty(loc::Real, scale::Real; interc::Real=0.0)

Struct representing a parametric penalty function.

# Fields
- `loc::Float64`: The location parameter.
- `scale::Float64`: The scale parameter.
- `interc::Float64`: The intercept parameter (defaults to 0.0).
"""
struct ParametricPenalty <: PenaltyFunction
    loc::Float64
    scale::Float64
    interc::Float64
    function ParametricPenalty(loc::Real, scale::Real; interc::Real=0.0)
        new(loc, scale, interc)
    end
end

"""
    (penalty::ParametricPenalty)(x::Vector{T}) where {T}

Applies the parametric penalty function to the input vector.

# Arguments
- `x::Vector{T}`: The input vector.

# Returns
- `Float64`: The calculated parametric penalty.
"""
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


