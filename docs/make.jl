using ONSAS
using Documenter

DocMeta.setdocmeta!(ONSAS, :DocTestSetup, :(using ONSAS); recursive=true)

makedocs(;
    modules=[ONSAS],
    authors="To be defined",
    repo="https://github.com/ONSAS/ONSAS.jl/blob/{commit}{path}#{line}",
    sitename="ONSAS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ONSAS.github.io/ONSAS.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Modules" => [
            "Materials" => [
                "AbstractMaterial" => "lib/Materials.md",
                "AbstractLinearElasticMaterial" => "lib/LinearElasticMaterials.md",
                "AbstractHyperElasticMaterials" => "lib/HyperElasticMaterials.md",
            ],
            "Elements" => "lib/Elements.md",
            "CrossSections" => "lib/CrossSections.md",
            "Meshes" => "lib/Meshes.md",
            "StructuralModel" => "lib/StructuralModel.md",
            "StructuralSolvers" => "lib/StructuralSolvers.md",
            "StructuralAnalysis" => "lib/StructuralAnalysis.md",
            "StaticAnalysis" => "lib/StaticAnalysis.md",
        ],
        "References" => "references.md",
    ]
)

deploydocs(;
    repo="github.com/ONSAS/ONSAS.jl",
    push_preview=true
)
