using NetworkDistances
using Documenter

DocMeta.setdocmeta!(NetworkDistances, :DocTestSetup, :(using NetworkDistances); recursive=true)

makedocs(
    modules = [
        NetworkDistances,
    ],
    authors = "George Bolt <g.bolt@lancaster.ac.uk> and contributors",
    repo = "https://github.com/your-username/NetworkDistances.jl/blob/{commit}{path}#{line}",
    sitename = "NetworkDistances.jl",
    format = Documenter.HTML(; 
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://your-username.github.io/NetworkDistances.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
    warnonly = [:missing_docs, :autodocs_block]
)

deploydocs(
    repo = "github.com/your-username/NetworkDistances.jl",
    devbranch = "main",
)


