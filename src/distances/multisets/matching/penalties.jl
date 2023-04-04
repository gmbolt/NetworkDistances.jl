export FixedPenalty, SizePenalty, DistancePenalty, ParametricPenalty

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
