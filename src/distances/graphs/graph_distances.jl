using LinearAlgebra, Distances
export Diffusion, Hamming
export hamming_dist, jaccard_dist, diffusion_dist

# Here we extend the Hamming and Jaccard distances defined in Distances.jl
# We also define a new distances: Diffusion. 

"""
    hamming_dist(A1::Matrix{Int}, A2::Matrix{Int})

Compute the Hamming distance between two adjacency matrices (directed graphs).

# Arguments
- `A1::Matrix{Int}`: The first adjacency matrix.
- `A2::Matrix{Int}`: The second adjacency matrix.

# Returns
- `Int`: The Hamming distance.
"""
function hamming_dist(A1::Matrix{Int}, A2::Matrix{Int})
    return sum(abs(x - y) for (x, y) in zip(A1, A2))
end

(d::Distances.Hamming)(A1::Matrix{Int}, A2::Matrix{Int}) = hamming_dist(A1, A2)

(d::Distances.Hamming)(A1::Matrix{Int}, A2::Nothing) = sum(A1)
(d::Distances.Hamming)(A1::Nothing, A2::Matrix{Int}) = d(A2, A1)


"""
    jaccard_dist(A1::Matrix{Int}, A2::Matrix{Int})

Compute the Jaccard distance between two adjacency matrices (directed graphs).

# Arguments
- `A1::Matrix{Int}`: The first adjacency matrix.
- `A2::Matrix{Int}`: The second adjacency matrix.

# Returns
- `Float64`: The Jaccard distance.
"""
function jaccard_dist(A1::Matrix{Int}, A2::Matrix{Int})
    return 1 - sum(min(x, y) for (x, y) in zip(A1, A2)) / sum(max(x, y) for (x, y) in zip(A1, A2))
end

(d::Distances.Jaccard)(A1::Matrix{Int}, A2::Matrix{Int}) = jaccard_dist(A1, A2)

"""
    Diffusion

Custom metric type for Diffusion distance.
"""
struct Diffusion <: Metric end

"""
    diffusion_dist(A1::Matrix{Int}, A2::Matrix{Int}, t)

Compute the diffusion distance between two adjacency matrices with `t` as the time of diffusion.

# Arguments
- `A1::Matrix{Int}`: The first adjacency matrix.
- `A2::Matrix{Int}`: The second adjacency matrix.
- `t`: The time of diffusion.

# Returns
- `Float64`: The diffusion distance.
"""
function diffusion_dist(A1::Matrix{Int}, A2::Matrix{Int}, t)
    L1 = Diagonal(vec(sum(A1, dims=2))) - A1
    L2 = Diagonal(vec(sum(A2, dims=2))) - A2
    function matrixExp(A, t)
        F = eigen(A)
        return F.vectors * Diagonal(exp.(F.values * t)) * transpose(F.vectors)
    end
    return sum((matrixExp(-L1, t) - matrixExp(-L2, t)) .^ 2)
end

(d::Diffusion)(A1::Matrix{Int}, A2::Matrix{Int}) = diffusion_dist(A1, A2)


