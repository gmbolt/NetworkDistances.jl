using RecipesBase
export PathSeqPlot, NoisePlot 
export get_entry_xlocs

"""
    get_entry_xlocs(lens::Vector{Int}, width, entrymargin, pathmargin; origin=0.0)

Calculates the x-coordinates for entries in a sequence of paths for plotting purposes.

# Arguments
- `lens::Vector{Int}`: A vector containing the lengths of individual paths.
- `width`: The width of each entry.
- `entrymargin`: The margin between entries within a path.
- `pathmargin`: The margin between different paths.
- `origin`: The starting x-coordinate (default: 0.0).

# Returns
- `Vector{Float64}`: A vector of x-coordinates for each entry.
"""
function get_entry_xlocs(
    lens::Vector{Int},
    width,
    entrymargin,
    pathmargin;
    origin=0.0
    )
    out = Float64[]
    x_tmp = origin
    for n in lens
        for i in 1:n 
            push!(out, x_tmp)
            x_tmp += width + entrymargin
        end 
        x_tmp -= entrymargin # Take of last entry margin
        x_tmp += pathmargin 
    end 
    return out
end 

"""
    PathSeqPlot(obs; entrymargin=0.1, entryfontsize=20, pathmargin=0.5, entrycolor=:green, align=:center)

Plots a sequence of paths, where each path is a sequence of observations.

# Arguments
- `obs`: The observations, typically a vector of vectors.
- `entrymargin`: Margin between entries within a path (default: 0.1).
- `entryfontsize`: Font size for annotations (default: 20).
- `pathmargin`: Margin between different paths (default: 0.5).
- `entrycolor`: Color of the entries (default: :green). Can be a single color or a vector of vectors for individual entry colors.
- `align`: Alignment of the plot (:center, :left, or :right, default: :center).
"""
@userplot PathSeqPlot
@recipe function f(
    h::PathSeqPlot; 
    entrymargin=0.1, 
    entryfontsize=20,
    pathmargin=0.5,
    entrycolor=:green,
    align=:center
    )
    obs = h.args[1]
    w,h = (1,1) 
    # Get entry centers 
    y = fill(h/2, sum(length, obs))
    x = get_entry_xlocs(
        length.(obs), 
        w, 
        entrymargin, 
        pathmargin, 
        origin=w/2
    )
    # Shift x according to alignment
    if align==:center 
        x .-= x[end]/2
    elseif align==:right 
        x .-= (x[end] + w/2)
    else
        if align!=:left 
            error("Align specification not supported. Must be, :center, :left or :right")
        end 
    end 
    annotations := [(x,y,string(i)) for (x,y,i) in zip(x,y,vcat(obs...))]
    annotationfontsize := entryfontsize
    x_cords, y_cords = rectangle_corners(x,y, w, h; anchor=:center)
    showaxis --> false 
    axis --> nothing 
    aspect_ratio --> 1 
    legend --> false 
    if typeof(entrycolor)==Symbol
        fillcolor --> entrycolor
    elseif typeof(entrycolor)<:Vector{Vector{T}} where {T}
        fillcolor --> permutedims(vcat(entrycolor...))
    end 
    @series begin
        seriestype := :shape
        x_cords,y_cords
    end 
end 

"""
    NoisePlot(obs, err; entrymargin=0.1, entryfontsize=20, stdcolor=:green, errcolor=:magenta, pathmargin=0.5, align=:center)

Plots a sequence of paths with an indication of errors or noise.

# Arguments
- `obs`: The observations, typically a vector of vectors.
- `err`: A boolean vector indicating error status for each observation.
- `entrymargin`: Margin between entries within a path (default: 0.1).
- `entryfontsize`: Font size for annotations (default: 20).
- `stdcolor`: Color for standard (non-error) entries (default: :green).
- `errcolor`: Color for error entries (default: :magenta).
- `pathmargin`: Margin between different paths (default: 0.5).
- `align`: Alignment of the plot (:center, :left, or :right, default: :center).
"""
@userplot NoisePlot
@recipe function f(
    h::NoisePlot; 
    entrymargin=0.1, 
    entryfontsize=20,
    stdcolor=:green, errcolor=:magenta,
    pathmargin=0.5,
    align=:center
    )
    obs, err = (h.args[1], h.args[2])
    w,h = (1,1) 
    # Get entry centers 
    y = fill(h/2, sum(length, obs))
    x = get_entry_xlocs(
        length.(obs), 
        w, 
        entrymargin, 
        pathmargin, 
        origin=w/2
    )
    # Shift x according to alignment
    if align==:center 
        x .-= x[end]/2
    elseif align==:right 
        x .-= (x[end] + w/2)
    else
        if align!=:left 
            error("Align specification not supported. Must be, :center, :left or :right")
        end 
    end 
    annotations := [(x,y,string(i)) for (x,y,i) in zip(x,y,vcat(obs...))]
    annotationfontsize := entryfontsize
    x_cords, y_cords = rectangle_corners(x,y, w, h; anchor=:center)
    showaxis --> false 
    axis --> nothing 
    aspect_ratio --> 1 
    legend --> false 
    fillcolor --> permutedims([i ? stdcolor : errcolor for i in vcat(err...)])
    @series begin
        seriestype := :shape
        x_cords,y_cords
    end 
end 


