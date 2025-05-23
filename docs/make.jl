using ONSAS
using Documenter

DocMeta.setdocmeta!(ONSAS, :DocTestSetup, :(using ONSAS); recursive = true)

makedocs(;
    modules = [ONSAS],
    repo = "https://github.com/ONSAS/ONSAS.jl/blob/{commit}{path}#{line}",
    sitename = "ONSAS.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://ONSAS.github.io/ONSAS.jl",
        edit_link = "main",
        assets = String[]),
    pages = ["Home" => "index.md",
        "Modules" => [
            "BoundaryConditions" => [
                "BoundaryConditions" => "lib/BoundaryConditions/BoundaryConditions.md",
                "DirichletBoundaryConditions" => "lib/BoundaryConditions/DirichletBoundaryConditions.md",
                "FixedFieldBoundaryConditions" => "lib/BoundaryConditions/FixedFieldBoundaryConditions.md",
                "GlobalLoadBoundaryConditions" => "lib/BoundaryConditions/GlobalLoadBoundaryConditions.md",
                "LocalLoadBoundaryConditions" => "lib/BoundaryConditions/LocalLoadBoundaryConditions.md"],
            "CrossSections" => ["Circles" => "lib/CrossSections/Circles.md",
                "CrossSections" => "lib/CrossSections/CrossSections.md",
                "GenericCrossSections" => "lib/CrossSections/GenericCrossSections.md",
                "Rectangles" => "lib/CrossSections/Rectangles.md",
                "Squares" => "lib/CrossSections/Squares.md"],
            "Entities" => ["Entities" => "lib/Entities/Entities.md",
                "Frames" => "lib/Entities/Frames.md",
                "Nodes" => "lib/Entities/Nodes.md",
                "Tetrahedrons" => "lib/Entities/Tetrahedrons.md",
                "TriangularFaces" => "lib/Entities/TriangularFaces.md",
                "Trusses" => "lib/Entities/Trusses.md"],
            "Interfaces" => ["Gmsh" => "lib/Interfaces/Gmsh.md",
                "VTK" => "lib/Interfaces/VTK.md"],
            "Materials" => [
                "HyperElasticMaterial" => "lib/Materials/HyperElasticMaterial.md",
                "HyperElasticMaterials" => "lib/Materials/HyperElasticMaterials.md",
                "IsotropicLinearElasticMaterial" => "lib/Materials/IsotropicLinearElasticMaterial.md",
                "LinearElasticMaterials" => "lib/Materials/LinearElasticMaterials.md",
                "Materials" => "lib/Materials/Materials.md",
                "NeoHookeanMaterial" => "lib/Materials/NeoHookeanMaterial.md",
                "SVKMaterial" => "lib/Materials/SVKMaterial.md"],
            "Meshes" => ["Handlers" => "lib/Meshes/Handlers.md",
                "Interpolators" => "lib/Meshes/Interpolators.md",
                "Meshes" => "lib/Meshes/Meshes.md",
                "Searches" => "lib/Meshes/Searches.md"],
            "StructuralAnalyses" => [
                "LinearStaticAnalyses" => "lib/StructuralAnalyses/LinearStaticAnalyses.md",
                "NonLinearStaticAnalyses" => "lib/StructuralAnalyses/NonLinearStaticAnalyses.md",
                "StaticAnalyses" => "lib/StructuralAnalyses/StaticAnalyses.md",
                "StaticStates" => "lib/StructuralAnalyses/StaticStates.md",
                "StructuralAnalyses" => "lib/StructuralAnalyses/StructuralAnalyses.md"],
            "StructuralModel" => [
                "StructuralBoundaryConditions" => "lib/StructuralModel/StructuralBoundaryConditions.md",
                "StructuralEntities" => "lib/StructuralModel/StructuralEntities.md",
                "StructuralMaterials" => "lib/StructuralModel/StructuralMaterials.md",
                "Structures" => "lib/StructuralModel/Structures.md"],
            "StructuralSolvers" => ["Assemblers" => "lib/StructuralSolvers/Assemblers.md",
                "Solutions" => "lib/StructuralSolvers/Solutions.md",
                "Solvers" => "lib/StructuralSolvers/Solvers.md",
                "StructuralSolvers" => "lib/StructuralSolvers/StructuralSolvers.md"],
            "Utils" => ["Utils" => "lib/Utils/Utils.md"]],
        "References" => "references.md"])

deploydocs(;
    repo = "github.com/ONSAS/ONSAS.jl",
    push_preview = true)
