using ONSAS
using Documenter

DocMeta.setdocmeta!(ONSAS, :DocTestSetup, :(using ONSAS); recursive=true)

makedocs(;
    modules=[ONSAS],
    authors="Jorge Pérez Zerpa",
    repo="https://github.com/ONSAS/ONSAS.jl/blob/{commit}{path}#{line}",
    sitename="ONSAS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ONSAS.github.io/ONSAS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ONSAS/ONSAS.jl",
    devbranch="main",
)
