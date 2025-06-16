using Hungarian

# Distance evaluation 
# -------------------

# Note each matching-based distance must have two methods 
# 1. get_cost_matrix_dynamic() - this will swap the dimensions of the cost matrix depending on size of passed objects
# 2. get_cost_matrix_fixed() - this will keep the demension fixed, so that the dimensions of cost matrix will depend on size of passed object

Base.show(io::IO, d::CompleteMatchingDistance) = print(io, typeof(d))

"""
    (d::T)(X::Vector{S}, Y::Vector{S}) where {T<:CompleteMatchingDistance, S}

Computes the complete matching distance between two vectors `X` and `Y`.

# Arguments
- `d::T`: The complete matching distance metric.
- `X::Vector{S}`: The first vector of elements.
- `Y::Vector{S}`: The second vector of elements.

# Returns
- `Float64`: The computed matching distance.
"""
function (d::T where {T<:CompleteMatchingDistance})(
    X::Vector{S}, Y::Vector{S}
) where {S}
    C = get_cost_matrix_dynamic(d, X, Y)
    return eval_distance(d.optimiser, C)
end

function Base.show(io::IO, d_gen::General{T}) where {T<:CompleteMatchingDistance}
    b = IOBuffer()
    show(b, d_gen.d)
    d_str = String(take!(b))
    print(io, "General{$(d_str)}")
end

"""
    (d::General{T})(X::Vector{S}, Y::Vector{S}) where {T<:CompleteMatchingDistance, S}

Computes the general matching distance between two vectors `X` and `Y`.

# Arguments
- `d::General{T}`: The general matching distance metric.
- `X::Vector{S}`: The first vector of elements.
- `Y::Vector{S}`: The second vector of elements.

# Returns
- `Float64`: The computed general matching distance.
"""
function (d::General{T} where {T<:CompleteMatchingDistance})(
    X::Vector{S}, Y::Vector{S}
) where {S}
    C = get_cost_matrix_fixed(d, X, Y)
    return hungarian(C)[2]
end


