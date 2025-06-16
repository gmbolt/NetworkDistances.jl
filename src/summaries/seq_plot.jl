using RecipesBase
export SeqPlot
export TestPlot
export rectangle_corners

"""
    rectangle_corners(x::Real,y::Real,w,h; anchor=:bottomright)

Calculates the corner coordinates for a single rectangle.

# Arguments
- `x::Real`: The x-coordinate of the anchor point.
- `y::Real`: The y-coordinate of the anchor point.
- `w`: The width of the rectangle.
- `h`: The height of the rectangle.
- `anchor::Symbol`: The anchor point of the rectangle (:bottomright or :center).

# Returns
- `Tuple{Vector{Real}, Vector{Real}}`: A tuple containing two vectors, one for x-coordinates and one for y-coordinates of the rectangle corners.
"""
function rectangle_corners(x::Real,y::Real,w,h; anchor=:bottomright)
    if anchor == :botttomright 
        [x,x+w,x+w,x], [y,y,y+h,y+h]
    elseif anchor == :center 
        [x-w/2, x+w/2, x+w/2, x-w/2], [y-h/2, y-h/2, y+h/2, y+h/2]
    else 
        error("Anchor not recognised.")
    end 
end 

"""
    rectangle_corners(x_vec::Vector{T},y_vec::Vector{T},w,h; anchor=:bottomleft) where {T<:Real}

Calculates the corner coordinates for multiple rectangles based on input vectors of anchor points.

# Arguments
- `x_vec::Vector{T}`: A vector of x-coordinates for the anchor points.
- `y_vec::Vector{T}`: A vector of y-coordinates for the anchor points.
- `w`: The width of each rectangle.
- `h`: The height of each rectangle.
- `anchor::Symbol`: The anchor point of the rectangles (:bottomleft or :center).

# Returns
- `Tuple{Vector{Vector{Real}}, Vector{Vector{Real}}}`: A tuple containing two vectors of vectors, one for x-coordinates and one for y-coordinates of the rectangle corners.
"""
function rectangle_corners(
    x_vec::Vector{T},y_vec::Vector{T},
    w,h; anchor=:bottomleft
    ) where {T<:Real}

    x_out = Vector{Real}[]
    y_out = Vector{Real}[]
    if anchor == :bottomleft
        for (x,y) in zip(x_vec, y_vec) 
            push!(x_out, [x,x+w,x+w,x])
            push!(y_out, [y,y,y+h,y+h])
        end
    elseif anchor == :center 
        for (x,y) in zip(x_vec, y_vec) 
            push!(x_out, [x-w/2, x+w/2, x+w/2, x-w/2])
            push!(y_out, [y-h/2, y-h/2, y+h/2, y+h/2])
        end
    else 
        error("Anchor not recognised.")
    end 
    x_out, y_out
end 

"""
    SeqPlot(seq; entrymargin=0.1, entryfontsize=20)

Plots a sequence of elements as a series of rectangles with annotations.

# Arguments
- `seq`: The sequence of elements to plot.
- `entrymargin`: Margin between entries (default: 0.1).
- `entryfontsize`: Font size for annotations (default: 20).
"""
@userplot SeqPlot
@recipe function f(
    h::SeqPlot; 
    entrymargin=0.1, 
    entryfontsize=20
    )
    seq = h.args[1]
    w,h = (1,1) 
    y = fill(0.0, length(seq))
    x = [0.0 + i*(w+entrymargin) for i in 0:(length(seq)-1)]
    annotations := [(x,y,string(i)) for (x,y,i) in zip(x,y,seq)]
    annotationfontsize := entryfontsize
    x_cords, y_cords = rectangle_corners(x,y, w, h; anchor=:center)
    showaxis --> false 
    axis --> nothing 
    aspect_ratio --> 1 
    legend --> false 
    @series begin
        seriestype := :shape
        x_cords,y_cords
    end 
end 

"""
    TestPlot()

A test plot recipe for demonstration purposes.
"""
@userplot TestPlot
@recipe function f(
    h::TestPlot
    )
    aspect_ratio --> 1
    showaxis --> false 
    @series begin 
        seriestype := :shape
        1:4,1:4
    end 
    @series begin
        seriestype := :line
        1:4, 1:4
        annotations := [[];[(i,i,string(i)) for i in 1:4]]
        annotationhalign --> :right
    end
end 


