"""
    emd(a::AbstractVector, b::AbstractVector, C::AbstractMatrix)

Wrapper for `ot.emd()` method of the POT (Python Optimal Transport) python package.
Computes the Earth Mover's Distance (EMD) and returns the optimal transportation plan.

# Arguments
- `a::AbstractVector`: The weights of the first distribution.
- `b::AbstractVector`: The weights of the second distribution.
- `C::AbstractMatrix`: The cost matrix between the two distributions.

# Returns
- `Matrix{Float64}`: The optimal transportation plan.
"""
function emd(a::AbstractVector, b::AbstractVector, C::AbstractMatrix)
    return pyconvert(
        Matrix{Float64},
        ot[].emd(
            np[].array(a), np[].array(b),
            np[].array(C)
        )
    )
end

"""
    emd2(a::AbstractVector, b::AbstractVector, C::AbstractMatrix)

Wrapper for `ot.emd2()` method of the POT (Python Optimal Transport) python package.
Computes the Earth Mover's Distance (EMD) and returns the cost.

# Arguments
- `a::AbstractVector`: The weights of the first distribution.
- `b::AbstractVector`: The weights of the second distribution.
- `C::AbstractMatrix`: The cost matrix between the two distributions.

# Returns
- `Float64`: The Earth Mover's Distance cost.
"""
function emd2(a::AbstractVector, b::AbstractVector, C::AbstractMatrix)
    return pyconvert(
        Float64,
        ot[].emd2(
            np[].array(a), np[].array(b),
            np[].array(C)
        )
    )
end


