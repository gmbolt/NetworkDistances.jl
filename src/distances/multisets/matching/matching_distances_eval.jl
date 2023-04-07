using PythonOT, Hungarian

# Distance evaluation 
# -------------------

# Note each matching-based distance must have two methods 
# 1. get_cost_matrix_dynamic() - this will swap the dimensions of the cost matrix depending on size of passed objects
# 2. get_cost_matrix_fixed() - this will keep the demension fixed, so that the dimensions of cost matrix will depend on size of passed object

Base.show(io::IO, d::CompleteMatchingDistance) = print(io, typeof(d))

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

function (d::General{T} where {T<:CompleteMatchingDistance})(
    X::Vector{S}, Y::Vector{S}
) where {S}
    C = get_cost_matrix_fixed(d, X, Y)
    return hungarian(C)[2]
end
