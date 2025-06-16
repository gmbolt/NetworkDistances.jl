using NetworkLayout, Roots
export MultigraphPlot, relabel, get_curves

"""
    relabel(edgelist::Vector{Tuple{Int,Int}})

Relabels integer-based edge lists to a new contiguous integer range starting from 1.
It also returns a reverse mapping to the original labels.

# Arguments
- `edgelist::Vector{Tuple{Int,Int}}`: A vector of tuples representing edges with integer labels.

# Returns
- `Tuple{Vector{Tuple{Int,Int}}, Dict{Int,Int}}`: A tuple containing the relabeled edge list and a dictionary for reverse mapping.
"""
function relabel(edgelist::Vector{Tuple{Int,Int}})
    edgelist_relab = Tuple{Int,Int}[]
    vmap = Dict{Int,Int}()
    i = 1
    for edge in edgelist
        for v in edge
            if v ∉ keys(vmap)
                vmap[v] = i 
                i+=1
            end 
        end 
        push!(edgelist_relab, map(x->vmap[x], edge))
    end 
    rev_vmap = Dict(v=>k for (k,v) in vmap)
    return edgelist_relab, rev_vmap
end 

"""
    relabel(edgelist::Vector{Tuple{String,String}})

Relabels string-based edge lists to a new contiguous integer range starting from 1.
It also returns a reverse mapping to the original labels.

# Arguments
- `edgelist::Vector{Tuple{String,String}}`: A vector of tuples representing edges with string labels.

# Returns
- `Tuple{Vector{Tuple{Int,Int}}, Dict{Int,String}}`: A tuple containing the relabeled edge list and a dictionary for reverse mapping.
"""
function relabel(edgelist::Vector{Tuple{String,String}})
    edgelist_relab = Tuple{Int,Int}[]
    vmap = Dict{String,Int}()
    i = 1
    for edge in edgelist
        for v in edge
            if v ∉ keys(vmap)
                vmap[v] = i 
                i+=1
            end 
        end 
        push!(edgelist_relab, map(x->vmap[x], edge))
    end 
    rev_vmap = Dict(v=>k for (k,v) in vmap)
    return edgelist_relab, rev_vmap
end 

"""
    node_edge_intercept(x,w,h,r)

Calculates the intercept for a node edge.
"""
node_edge_intercept(x,w,h,r) = (4*h/(w^2) * x*(w-x))^2 + x^2 - r^2

"""
    get_arc(src, dst; bendprop=0.1)

Generates coordinates for a curved arc between two points.

# Arguments
- `src`: Source point coordinates.
- `dst`: Destination point coordinates.
- `bendprop`: Proportion of bending for the arc (default: 0.1).

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: X and Y coordinates of the arc.
"""
function get_arc(
    src, dst;
    bendprop=0.1
    )
    # f(x) = 4h/w² x(w-x) (equation for arc from (0,0)->(w,0))
    # We then then translate and rotate this 
    w = euclidean(src,dst)
    h = bendprop*euclidean(src,dst)
    θ = angle(complex(dst...)-complex(src...))
    cosθ, sinθ = (cos(θ), sin(θ))
    C = 4*h/w^2
    rawx = 0:0.01:w
    srcx,srcy = src
    xvals = [x*cosθ - C*x*(w-x)*sinθ + srcx for x in rawx]
    yvals = [x*sinθ + C*x*(w-x)*cosθ + srcy for x in rawx]
    return xvals, yvals
end

"""
    get_arc_shorten(src, dst; bendprop=0.1, shorten=0.2)

Generates coordinates for a shortened curved arc between two points.

# Arguments
- `src`: Source point coordinates.
- `dst`: Destination point coordinates.
- `bendprop`: Proportion of bending for the arc (default: 0.1).
- `shorten`: Amount to shorten the arc from both ends (default: 0.2).

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: X and Y coordinates of the shortened arc.
"""
function get_arc_shorten(
    src, dst;
    bendprop=0.1,
    shorten=0.2
    )

    w = euclidean(src,dst)
    h = bendprop*euclidean(src,dst)
    x_shorten = find_zero(x->node_edge_intercept(x,w,h,shorten), (0.0,shorten))
    θ = angle(complex(dst...)-complex(src...))
    cosθ, sinθ = (cos(θ), sin(θ))
    C = 4*h/w^2
    rawx = x_shorten:0.01:(w-x_shorten)
    srcx,srcy = src
    xvals = [x*cosθ - C*x*(w-x)*sinθ + srcx for x in rawx]
    yvals = [x*sinθ + C*x*(w-x)*cosθ + srcy for x in rawx]
    return xvals, yvals
end

"""
    get_curves(edges::AbstractVector, locs::AbstractVector; bendprop=0.1, shorten=0.2)

Generates curves for a list of edges given node locations.

# Arguments
- `edges::AbstractVector`: A list of edges (tuples of node indices).
- `locs::AbstractVector`: A list of node coordinates.
- `bendprop`: Proportion of bending for the arcs (default: 0.1).
- `shorten`: Amount to shorten the arcs from both ends (default: 0.2).

# Returns
- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Float64}}}`: X and Y coordinates for all curves.
"""
function get_curves(
    edges::AbstractVector,
    locs::AbstractVector;
    bendprop=0.1,
    shorten=0.2
    )

    xout, yout = (Vector{Float64}[],Vector{Float64}[])

    for (i,j) in edges
        xtmp, ytmp = get_arc_shorten(locs[i],locs[j],bendprop=bendprop,shorten=shorten)
        push!(xout,xtmp)
        push!(yout,ytmp)
    end 

    return xout, yout

end 

"""
    MultigraphPlot(edgelist; relabelnodes=true, shorten=0.1, bendprop=0.1, seed=1)

Plots a multigraph with curved edges.

# Arguments
- `edgelist`: A list of edges, either as tuples of integers or strings.
- `relabelnodes`: Whether to relabel nodes to contiguous integers (default: `true`).
- `shorten`: Amount to shorten the edges from nodes (default: 0.1).
- `bendprop`: Proportion of bending for the edges (default: 0.1).
- `seed`: Seed for the graph layout algorithm (default: 1).
"""
@userplot MultigraphPlot
@recipe function f(
    h::MultigraphPlot; 
    relabelnodes=true,
    shorten=0.1,
    bendprop=0.1,
    seed=1
    )
    edgelist = h.args[1]
    
    do_relabel = (relabelnodes) || (eltype(edgelist)==Tuple{String,String})
    if do_relabel
        edgelist, rev_vmap = relabel(edgelist)
        V = length(keys(rev_vmap))
        adjmat = zeros(Int, V, V)
        for (i,j) in edgelist
            adjmat[i,j] += 1 
        end 
    end 
    aspect_ratio --> 1
    locs = NetworkLayout.spring(adjmat,seed=seed)
    
    @series begin 
        # seriestype := :line
        linecolor := :blue
        arrow --> true
        get_curves(edgelist, locs, shorten=shorten, bendprop=bendprop)
    end 

    @series begin 
        seriestype := :scatter
        markercolor := :blue
        locs
    end 
end 

