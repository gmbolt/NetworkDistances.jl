# NetworkDistances.jl

NetworkDistances.jl is a Julia package for computing various distances between networks, paths, multisets, and sequences. It leverages `PythonCall.jl` for Earth Mover's Distance calculations via the Python Optimal Transport (POT) library.

## Installation

To install NetworkDistances.jl, you need to have Julia installed. Then, you can add the package from the Julia package manager:

```julia
using Pkg
Pkg.add("NetworkDistances")
```

## Usage

### Graph Distances

```julia
using NetworkDistances, Distances

A1 = [0 1; 1 0] # Adjacency matrix for a simple graph
A2 = [0 1; 0 0] # Another adjacency matrix

# Hamming distance
hamming_dist(A1, A2)

# Jaccard distance
jaccard_dist(A1, A2)

# Diffusion distance
diffusion_dist(A1, A2, 0.5) # with diffusion time t = 0.5
```

### Multiset Distances (Earth Mover's Distance)

```julia
using NetworkDistances, Distances

X = [1, 2, 3]
Y = [2, 3, 4]

# Earth Mover's Distance with AbsoluteDiff as ground distance
emd_dist = EMD(AbsoluteDiff())
emd_dist(X, Y)

# Scaled Earth Mover's Distance
semd_dist = sEMD(Euclidean(), AbsoluteDiff(), 0.5)
semd_dist(X, Y)
```

### Path Distances (LCS and LSP)

```julia
using NetworkDistances

X = [1, 2, 3, 4]
Y = [1, 3, 2, 4]

# Longest Common Subsequence (LCS) distance
lcs_dist = LCS()
lcs_dist(X, Y)

# Longest Common Subpath (LSP) distance
lsp_dist = LSP()
lsp_dist(X, Y)
```

### Pre-computed Distances

```julia
using NetworkDistances, Distances

data = [1, 2, 3, 4]
precomputed_euclidean = pre_compute(Euclidean(), data)
precomputed_euclidean(1, 2)
```

### Set Distances

```julia
using NetworkDistances, Distances

set1 = Set([1, 2, 3])
set2 = Set([2, 3, 4])

# Hamming distance for sets
Distances.Hamming()(set1, set2)

# Jaccard distance for sets
Distances.Jaccard()(set1, set2)
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


