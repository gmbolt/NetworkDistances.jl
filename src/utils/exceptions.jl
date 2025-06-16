# Defines custom exception types for throwing errors 
"""
    NotImplementedError

Custom exception type for indicating that a feature or method is not yet implemented.

# Fields
- `msg::String`: A message describing the unimplemented feature.
"""
struct NotImplementedError <: Exception
    msg::String
end


