using StatsBase, Distances

export EarthMoversDistance, EMD
export check_trans_plan, get_info

# Optimal Transport (OT) Distances 
# --------------------------------

"""
    LengthDistance

Abstract type for defining distance metrics between lengths.
"""
abstract type LengthDistance <: SemiMetric end

"""
    AbsoluteDiff

Concrete type for `LengthDistance` representing the absolute difference between two integers.
"""
struct AbsoluteDiff <: LengthDistance end
function (d::AbsoluteDiff)(N::Int, M::Int)
    return abs(N - M)
end

"""
    SquaredDiff

Concrete type for `LengthDistance` representing the squared difference between two integers.
"""
struct SquaredDiff <: LengthDistance end
function (d::SquaredDiff)(N::Int, M::Int)
    return (N - M)^2
end

"""
    EarthMoversDistance{T<:SemiMetric}

Struct representing the Earth Mover's Distance (EMD).

# Fields
- `ground_dist::T`: The ground distance metric used for EMD.
"""
struct EarthMoversDistance{T<:SemiMetric} <: SemiMetric
    ground_dist::T
end

const EMD = EarthMoversDistance

"""
    check_trans_plan(d::EMD, X::Vector{T}, Y::Vector{T}) where {T}

Checks the validity of a transportation plan for Earth Mover's Distance.

# Arguments
- `d::EMD`: The Earth Mover's Distance metric.
- `X::Vector{T}`: The first set of data.
- `Y::Vector{T}`: The second set of data.

# Returns
- `Matrix{Float64}`: The cost matrix `C` used in the EMD calculation.

# Throws
- `error`: If the sum of the transportation plan is not approximately 1.0.
"""
function check_trans_plan(d::EMD, X::Vector{T}, Y::Vector{T}) where {T}

    a = proportionmap(X)
    b = proportionmap(Y)

    C = zeros(Float64, length(a), length(b))

    for (j, val_b) in enumerate(keys(b))
        for (i, val_a) in enumerate(keys(a))
            C[i, j] = d.ground_dist(val_a, val_b)
        end
    end

    TransPlan::Array{Float64,2} = emd(
        collect(values(a)),
        collect(values(b)), C
    )

    if sum(TransPlan) ≈ 1.0
        println("Transplan looks valid!")
    else
        error("Transplan looks spurious!")
    end
    return C
end

"""
    (d::EMD)(X::Vector{T}, Y::Vector{T}) where {T}

Computes the Earth Mover's Distance between two vectors.

# Arguments
- `d::EMD`: The Earth Mover's Distance metric.
- `X::Vector{T}`: The first vector.
- `Y::Vector{T}`: The second vector.

# Returns
- `Float64`: The Earth Mover's Distance.
"""
function (d::EMD)(X::Vector{T}, Y::Vector{T}) where {T}

    a = proportionmap(X)
    b = proportionmap(Y)

    C = zeros(Float64, length(a), length(b))

    for (j, val_b) in enumerate(keys(b))
        for (i, val_a) in enumerate(keys(a))
            C[i, j] = d.ground_dist(val_a, val_b)
        end
    end
    # Function emd2() defined in NetworkDistances.jl - uses calls of 
    # python functions via PythonCall.jl
    return emd2(
        collect(values(a)),
        collect(values(b)), C
    )
end

(d::EMD)(X::Vector{T}, Y::Nothing) where {T} = mean(d.ground_dist(xi, nothing) for xi in X)
(d::EMD)(X::Nothing, Y::Vector{T}) where {T} = d(Y, X)


"""
    get_info(d::EMD, X::Vector{T}, Y::Vector{T}) where {T}

Retrieves information about the Earth Mover's Distance calculation, including the keys of the proportion maps and the transportation plan.

# Arguments
- `d::EMD`: The Earth Mover's Distance metric.
- `X::Vector{T}`: The first vector.
- `Y::Vector{T}`: The second vector.

# Returns
- `Tuple`: A tuple containing the keys of proportion map `a`, keys of proportion map `b`, and the transportation plan.
"""
function get_info(
    d::EMD,
    X::Vector{T}, Y::Vector{T}
) where {T}

    a = proportionmap(X)
    b = proportionmap(Y)

    C = zeros(Float64, length(a), length(b))

    for (j, val_b) in enumerate(keys(b))
        for (i, val_a) in enumerate(keys(a))
            C[i, j] = d.ground_dist(val_a, val_b)
        end
    end

    TransPlan = emd(
        collect(values(a)),
        collect(values(b)), C
    )

    return collect(keys(a)), collect(keys(b)), TransPlan
end


"""
    sEMD{T<:SemiMetric,G<:LengthDistance}

Struct representing a scaled Earth Mover's Distance (sEMD) composed with a length distance.

# Fields
- `ground_dist::T`: The ground distance metric for EMD.
- `length_dist::G`: The length distance metric.
- `τ::Real`: Relative weighting term (a proportion weighting EMD vs length distance, high τ => high EMD weighting).
"""
struct sEMD{T<:SemiMetric,G<:LengthDistance} <: SemiMetric
    ground_dist::T
    length_dist::G
    τ::Real # Relative weighting term (a proportion weighting EMD vs length distance, high τ => high EMD weighting)
end

"""
    (d::sEMD)(S1::Vector{T}, S2::Vector{T}) where {T}

Computes the scaled Earth Mover's Distance between two vectors.

# Arguments
- `d::sEMD`: The scaled Earth Mover's Distance metric.
- `S1::Vector{T}`: The first vector.
- `S2::Vector{T}`: The second vector.

# Returns
- `Float64`: The scaled Earth Mover's Distance.
"""
function (d::sEMD)(S1::Vector{T}, S2::Vector{T}) where {T}

    d₁ = EMD(d.ground_dist)(S1, S2)
    d₂ = d.length_dist(length(S1), length(S2))
    # @show d₁, d₂

    return d₁ + d.τ * d₂


end

"""
    sEMD2{T<:SemiMetric,G<:LengthDistance}

Struct representing another variant of scaled Earth Mover's Distance (sEMD2).

# Fields
- `ground_dist::T`: The ground distance metric for EMD.
- `length_dist::G`: The length distance metric.
- `τ::Real`: Weighting term.
"""
struct sEMD2{T<:SemiMetric,G<:LengthDistance} <: SemiMetric
    ground_dist::T
    length_dist::G
    τ::Real
end

"""
    (d::sEMD2)(S1::Vector{T}, S2::Vector{T}) where {T}

Computes the second variant of scaled Earth Mover's Distance between two vectors.

# Arguments
- `d::sEMD2`: The second scaled Earth Mover's Distance metric.
- `S1::Vector{T}`: The first vector.
- `S2::Vector{T}`: The second vector.

# Returns
- `Float64`: The second scaled Earth Mover's Distance.
"""
function (d::sEMD2)(S1::Vector{T}, S2::Vector{T}) where {T}

    d₁ = EMD(d.ground_dist)(S1, S2)
    d₂ = d.length_dist(sum(length.(S1)), sum(length.(S2)))
    # @show d₁, d₂

    return d.τ * d₁ + (1 - d.τ) * d₂


end


