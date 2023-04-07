export get_info, print_info, get_info_deep

function deep_invert(x::Vector{Vector{Bool}})
    return [.!xi for xi in x]
end


# Matching distances 
# ==================

const MatchingBasedDistance = Union{CompleteMatchingDistance,General{T}} where {T<:CompleteMatchingDistance}

function print_info(
    d::T,
    X::Vector{S}, Y::Vector{S}
) where {T<:MatchingBasedDistance,S}

    title = "Matching Summary"
    println(title)
    println("-"^length(title) * "\n")
    indx, indy = get_info(d, X, Y)
    max_len = maximum(map(z -> length(@sprintf("%s", z)), X[indx]))
    for (i, j) in zip(indx, indy)
        xtmp, ytmp = (X[i], Y[j])
        pad = max_len - length(@sprintf("%s", xtmp))
        println(" "^pad * "$xtmp", " → $ytmp")
    end
    println("\nUnmatched entries (first):")
    for i in eachindex(X)
        if i ∉ indx
            print("$(X[i]), ")
        end
    end
    println("\nUnmatched entries (second):")
    for i in eachindex(Y)
        if i ∉ indy
            print("$(Y[i]), ")
        end
    end
end


function get_info(
    d::T,
    X::Vector{S}, Y::Vector{S}
) where {T<:MatchingBasedDistance,S}

    C = get_cost_matrix_fixed(d, X, Y)
    assignment, cost = hungarian(C)
    indx, indy = (Int[], Int[])
    for i in eachindex(X)
        j = assignment[i]
        if j ≤ length(Y)
            push!(indx, i)
            push!(indy, j)
        end
    end
    return indx, indy
end

"""
    get_info_deep(d::MatchingBasedDistance, X::Vector{Vector{S}}, Y::Vector{Vector{S}}) 

Returns two vectors indicating which entries have been matched. An entry being `true` in the output implies this entry was matched with one in the other observation. 
"""
function get_info_deep(
    d::T,
    X::Vector{S}, Y::Vector{S};
    invert::Bool=false
)::Tuple{Vector{BitVector},Vector{BitVector}} where {T<:MatchingBasedDistance,S}

    d_g = d.ground_dist

    indx, indy = get_info(d, X, Y)

    outx, outy = (
        [zeros(Bool, i) for i in length.(X)],
        [zeros(Bool, i) for i in length.(Y)]
    )

    for (i, j) in zip(indx, indy)
        outx[i], outy[j] = get_info(d_g, X[i], Y[j])
    end

    if invert
        return deep_invert(outx), deep_invert(outy)
    else
        return outx, outy
    end
end

# Edit distances
# ==============



function print_info(
    d::Union{EditDistance,FastEditDistance,FixPenEditDist,FastFixPenEditDist},
    x::Vector{T}, y::Vector{T}
) where {T}

    indx, indy = get_info(d, x, y)
    max_len = maximum(map(z -> length(@sprintf("%s", z)), x))
    ix, iy = (0, 0)
    title = "\nOptimal Matching"
    println(title)
    println("-"^length(title), "\n")
    while true
        ixerr = isnothing(findnext(indx, ix + 1)) ? length(x) : findnext(indx, ix + 1) - 1
        iyerr = isnothing(findnext(indy, iy + 1)) ? length(y) : findnext(indy, iy + 1) - 1
        for j in (ix+1):ixerr
            tmp1 = x[j]
            tmp2 = "Null"
            pad = max_len - length(@sprintf("%s", tmp1))
            println(" "^pad * "$tmp1", " → $tmp2")
        end
        for j in (iy+1):iyerr
            tmp1 = "Null"
            tmp2 = y[j]
            pad = max_len - length(@sprintf("%s", tmp1))
            println(" "^pad * "$tmp1", " → $tmp2")
        end
        if isnothing(findnext(indx, ix + 1))
            break
        end
        ix = findnext(indx, ix + 1)
        iy = findnext(indy, iy + 1)
        tmp1, tmp2 = (x[ix], y[iy])
        pad = max_len - length(@sprintf("%s", tmp1))
        println(" "^pad * "$tmp1", " → $tmp2")
    end
end



function get_info(
    d::Union{FixPenEditDist,FastFixPenEditDist},
    x::Vector{T}, y::Vector{T}
) where {T}
    C = zeros(Float64, length(x) + 1, length(y) + 1)

    C[:, 1] = [d.ρ * i for i = 0:length(x)]
    C[1, :] = [d.ρ * i for i = 0:length(y)]

    for j = 1:length(y)
        for i = 1:length(x)
            C[i+1, j+1] = minimum([
                C[i, j] + d.ground_dist(x[i], y[j]),
                C[i, j+1] + d.ρ,
                C[i+1, j] + d.ρ
            ])
        end
    end
    # Now retrace steps to determine an optimal matching
    i, j = size(C)
    indx, indy = (zeros(Bool, length(x)), zeros(Bool, length(y)))
    while (i ≠ 1) | (j ≠ 1)
        if C[i, j] == (C[i-1, j] + d.ρ)
            i = i - 1
        elseif C[i, j] == (C[i, j-1] + d.ρ)
            j = j - 1
        else
            i = i - 1
            j = j - 1
            indx[i], indy[j] = (true, true)
        end
    end
    return indx, indy
end

function get_info_deep(
    d::Union{EditDistance,FastEditDistance,FixPenEditDist,FastFixPenEditDist},
    x::Vector{T}, y::Vector{T};
    invert::Bool=false
)::Tuple{Vector{BitVector},Vector{BitVector}} where {T}

    d_g = d.ground_dist
    indx, indy = get_info(d, x, y)
    outx, outy = (
        [zeros(Bool, i) for i in length.(x)],
        [zeros(Bool, i) for i in length.(y)]
    )
    ix, iy = (0, 0)
    while true
        if isnothing(findnext(indx, ix + 1))
            break
        end
        ix = findnext(indx, ix + 1)
        iy = findnext(indy, iy + 1)
        outx[ix], outy[iy] = get_info(d_g, x[ix], y[iy])
    end

    if invert
        return deep_invert(outx), deep_invert(outy)
    else
        return outx, outy
    end
end