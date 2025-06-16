module NetworkDistances

using PythonCall
using LinearAlgebra, Distances, StatsBase

# References for python modules 
const np = Ref{Py}()    # Numpy
const ot = Ref{Py}()    # Python optimal transport

function __init__()
    np[] = pyimport("numpy")
    ot[] = pyimport("ot")
end

# Python Optimal Transport (POT) wrappers
include("py_wrappers.jl")

# Utilities 
include("utils/utils.jl")
include("utils/exceptions.jl")
include("utils/pairwise.jl")

# Path distances 
include("distances/paths/paths.jl")

# Graph distances 
include("distances/graphs/graph_distances.jl")

# Multiset distances
include("distances/multisets/matching/matching_distances_complete.jl")
include("distances/multisets/matching/matching_distances_generalised.jl")
include("distances/multisets/matching/matching_distances_eval.jl")
include("distances/multisets/earth_movers_distance.jl")

# Sequence distances
include("distances/sequences/edit_distance.jl")
include("distances/sequences/dtw_distance.jl")

# Pre-computed distances
include("distances/pre_computed.jl")

# Basic set distances (Hamming + Jaccard)
include("distances/set_distances.jl")

# Normalised distances (normalises any metric)
include("distances/normalised.jl")

# Summaries 
include("summaries/seq_plot.jl")
include("summaries/pathseq_plot.jl")
include("summaries/get_info.jl")
include("summaries/multigraph_plot.jl")

end


