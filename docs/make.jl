using NetworkDistances
using Documenter

DocMeta.setdocmeta!(NetworkDistances, :DocTestSetup, :(using NetworkDistances); recursive=true)

makedocs(;
    modules=[NetworkDistances],
    authors="George Bolt g.bolt@lancaster.ac.uk",
    repo="https://github.com/gmbolt/NetworkDistances.jl/blob/{commit}{path}#{line}",
    sitename="NetworkDistances.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gmbolt.github.io/NetworkDistances.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/gmbolt/NetworkDistances.jl"
    # devbranch="main"
)
