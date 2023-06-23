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
                                assets=String[]),
         pages=["Home" => "index.md",
                "Modules" => ["Materials" => ["AbstractMaterial" => "lib/Materials.md",
                                              "AbstractLinearElasticMaterial" => "lib/LinearElasticMaterials.md",
                                              "AbstractHyperElasticMaterials" => "lib/HyperElasticMaterials.md"],
                              "Entities" => "lib/Entities.md",
                              "CrossSections" => "lib/CrossSections.md",
                              "Meshes" => "lib/Meshes.md",
                              "Structures" => "lib/Structures.md",
                              "StructuralSolvers" => "lib/StructuralSolvers.md",
                              "StructuralAnalyses" => ["AbstractStructuralAnalyses" => "lib/StructuralAnalyses.md",
                                                       "AbstractStaticAnalysis" => ["AbstractStaticAnalysis" => "lib/StaticAnalyses.md",
                                                                                    "LinearStaticAnalysis" => "lib/LinearStaticAnalysis.md",
                                                                                    "NonLinearStaticAnalysis" => "lib/NonLinearStaticAnalysis.md"]]],
                "References" => "references.md"])

deploydocs(;
           repo="github.com/ONSAS/ONSAS.jl",
           push_preview=true)
