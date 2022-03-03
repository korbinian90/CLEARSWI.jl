using CLEARSWI
using Documenter

DocMeta.setdocmeta!(CLEARSWI, :DocTestSetup, :(using CLEARSWI); recursive=true)

makedocs(;
    modules=[CLEARSWI],
    authors="Korbinian Eckstein korbinian90@gmail.com",
    repo="https://github.com/korbinian90/CLEARSWI.jl/blob/{commit}{path}#{line}",
    sitename="CLEARSWI.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://korbinian90.github.io/CLEARSWI.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/korbinian90/CLEARSWI.jl",
    devbranch="master",
)
