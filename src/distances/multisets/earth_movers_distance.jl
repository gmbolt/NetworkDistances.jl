using StatsBase, Distances

export EarthMoversDistance, EMD
export check_trans_plan, get_info

# Optimal Transport (OT) Distances 
# --------------------------------

abstract type LengthDistance <: SemiMetric end

struct AbsoluteDiff <: LengthDistance end
function (d::AbsoluteDiff)(N::Int, M::Int)
    return abs(N - M)
end
struct SquaredDiff <: LengthDistance end
function (d::SquaredDiff)(N::Int, M::Int)
    return (N - M)^2
end

struct EarthMoversDistance{T<:SemiMetric} <: SemiMetric
    ground_dist::T
end

const EMD = EarthMoversDistance

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


# EMD composed with a distance of the number of interactions
struct sEMD{T<:SemiMetric,G<:LengthDistance} <: SemiMetric
    ground_dist::T
    length_dist::G
    τ::Real # Relative weighting term (a proportion weighting EMD vs length distance, high τ => high EMD weighting)
end

function (d::sEMD)(S1::Vector{T}, S2::Vector{T}) where {T}

    d₁ = EMD(d.ground_dist)(S1, S2)
    d₂ = d.length_dist(length(S1), length(S2))
    # @show d₁, d₂

    return d₁ + d.τ * d₂


end

struct sEMD2{T<:SemiMetric,G<:LengthDistance} <: SemiMetric
    ground_dist::T
    length_dist::G
    τ::Real
end

function (d::sEMD2)(S1::Vector{T}, S2::Vector{T}) where {T}

    d₁ = EMD(d.ground_dist)(S1, S2)
    d₂ = d.length_dist(sum(length.(S1)), sum(length.(S2)))
    # @show d₁, d₂

    return d.τ * d₁ + (1 - d.τ) * d₂


end

