module NetworkDistances
using PythonCall

# References for python modules 

# We use PythonCall.jl to allow use of the Python Optimal Transport (POT) (https://pythonot.github.io/)
# library for evaluating the earth mover's distance. Doing so within a Juli module/package
# requires first loading the necessary python modules/packages. This is done 
# as follows (following recommendations made in documentation for PythonCall.jl)...

# First define some references 
const np = Ref{Py}()    # Numpy
const ot = Ref{Py}()    # Python optimal transport

# Now using __init__() we can load required python modules into refs when 
# this Julia module is first loaded, as per recommendations of PythonCall.jl...
function __init__()
    np[] = pyimport("numpy")
    ot[] = pyimport("ot")
end

# With these modules imported and assigned to references, we can define some 
# wrappers to call the desired POT functions....
"""
Wrapper for ot.emd() method of POT python package.
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
Wrapper for ot.emd2() method of POT python package.
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

end
