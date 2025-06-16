using Distances

export Normalised
"""
    Normalised{T<:SemiMetric}

Normalised distance obtained by applying Steinhaus transform to a given distance `d`.

Note that this will be a metric if `d` is a metric.

# Fields
- `d::T`: The semi-metric to be normalised.
"""
struct Normalised{T<:SemiMetric} <: SemiMetric
    d::T
end


function Base.show(io::IO, d_n::Normalised{T}) where {T<:SemiMetric}
    b = IOBuffer()
    show(b, d_n.d)
    d_str = String(take!(b))
    print(io, "Normalised{$(d_str)}")
end

"""
    (d_n::Normalised{T})(x, y) where {T<:SemiMetric}

Computes the normalised distance between `x` and `y` using the Steinhaus transform.

# Arguments
- `d_n::Normalised{T}`: The normalised distance metric.
- `x`: The first element.
- `y`: The second element.

# Returns
- `Float64`: The normalised distance.
"""
function (d_n::Normalised{T})(x, y) where {T<:SemiMetric}
    d = d_n.d
    d_tmp = d(x, y)
    return 2 * d_tmp / (d(x, nothing) + d(y, nothing) + d_tmp)
end

(d_n::Normalised{T} where {T<:SemiMetric})(x::Nothing, y::Nothing) = 0.0 # To avoid NaN


