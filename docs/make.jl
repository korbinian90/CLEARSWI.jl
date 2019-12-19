using Documenter, SWI

makedocs(;
    modules=[SWI],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/korbinian90/SWI.jl/blob/{commit}{path}#L{line}",
    sitename="SWI.jl",
    authors="Korbinian Eckstein",
    assets=String[],
)

deploydocs(;
    repo="github.com/korbinian90/SWI.jl",
)
