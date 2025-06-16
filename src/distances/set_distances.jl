using Distances

"""
    (d::Distances.Hamming)(X::Set{T}, Y::Set{T}) where {T}

Computes the Hamming distance between two sets.

# Arguments
- `X::Set{T}`: The first set.
- `Y::Set{T}`: The second set.

# Returns
- `Int`: The Hamming distance.
"""
function (d::Distances.Hamming)(X::Set{T}, Y::Set{T}) where {T}
    if length(X) > length(Y) 
        return d(Y,X)
    else 
        n_intersect = 0
        for v in X 
            n_intersect += (v ∈ Y)
        end 
        return length(X) + length(Y) - 2n_intersect
    end 
end 

"""
    (d::Distances.Hamming)(X::Union{Set{T},Vector{T}}, Y::Nothing) where {T}
    (d::Distances.Hamming)(X::Nothing, Y::Union{Set{T},Vector{T}}) where {T}
    (d::Distances.Hamming)(X::Nothing,Y::Nothing)

Computes the Hamming distance when one or both inputs are `nothing`.
"""
(d::Distances.Hamming)(X::Union{Set{T},Vector{T}}, Y::Nothing) where {T} = length(X)
(d::Distances.Hamming)(X::Nothing, Y::Union{Set{T},Vector{T}}) where {T} = d(Y,X)
(d::Distances.Hamming)(X::Nothing,Y::Nothing) = 0.0

"""
    (d::Distances.Jaccard)(X::Set{T}, Y::Set{T}) where {T}

Computes the Jaccard distance between two sets.

# Arguments
- `X::Set{T}`: The first set.
- `Y::Set{T}`: The second set.

# Returns
- `Float64`: The Jaccard distance.
"""
function (d::Distances.Jaccard)(X::Set{T}, Y::Set{T}) where {T}
    if length(X) > length(Y) 
        return d(Y,X)
    else 
        n_intersect = 0
        for v in X 
            n_intersect += (v ∈ Y)
        end 

        return 1 - n_intersect/(length(X) + length(Y) - n_intersect)
    end 
end 


