# Matching distances

## Type hierarchy 

At the top level we have 

```julia
AbstractMatchingDistance <: SemiMetric
```
This encompasses all matching-based distances. Within this, we have two further abstract classes 
1. `CompleteMatchingDistance <: AbstractMatchingDistance` - these are the complete matching distances which optimise over *complete* matchings only (i.e. where all path from smaller observation are included in matching);
2. `General{T} <: AbstractMatchingDistance` where `T <: CompleteMatchingDistance` - these can have paths from both observations being unmatched. These are defined given complete matchings via `General(d)` where `d::CompleteMatchingDistance`.

## Other methods

```@docs
get_cost_matrix_dynamic
get_cost_matrix_fixed
```