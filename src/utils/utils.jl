using StatsBase
export get_aggregate_adj_mat

"""
    StatsBase.counts(data::Vector{T}, ref::Vector{T}) where {T}

Counts the occurrences of elements in `data` based on a reference `ref`.

# Arguments
- `data::Vector{T}`: The data vector.
- `ref::Vector{T}`: The reference vector containing unique elements.

# Returns
- `Vector{Int}`: A vector where each element is the count of the corresponding element in `ref` found in `data`.
"""
function StatsBase.counts(data::Vector{T}, ref::Vector{T}) where {T}
    mapper = Dict(val=>i for (i,val) in enumerate(ref))
    out = zeros(Int, length(ref))
    for x in data 
        out[mapper[x]] += 1
    end 
    return out
end 

"""
    get_aggregate_adj_mat(data::Vector{Vector{T}}, ref::Vector{T}) where {T}

Computes an aggregate adjacency matrix from a collection of paths.

# Arguments
- `data::Vector{Vector{T}}`: A vector of paths, where each path is a vector of elements.
- `ref::Vector{T}`: The reference vector containing unique elements that define the nodes of the adjacency matrix.

# Returns
- `Matrix{Int}`: An aggregate adjacency matrix where `A[i, j]` represents the number of times an edge from `ref[i]` to `ref[j]` appears in the paths.
"""
function get_aggregate_adj_mat(data::Vector{Vector{T}}, ref::Vector{T}) where {T}
    mapper = Dict(val=>i for (i,val) in enumerate(ref))
    A = zeros(Int, length(ref), length(ref))
    for path in data
        for i in Iterators.rest(eachindex(path),1)
            A[mapper[path[i-1]], mapper[path[i]]] += 1
        end 
    end 
    return A
end 


