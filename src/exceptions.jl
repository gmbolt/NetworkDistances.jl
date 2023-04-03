# Defines custom exception types for throwing errors 
struct NotImplementedError <: Exception
    msg::String
end
