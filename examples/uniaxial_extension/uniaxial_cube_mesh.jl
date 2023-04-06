using Gmsh

"Creates a mesh for a cube with a loaded face and a fixed origin at (0,0,0)
and (Lᵢ,Lⱼ,Lₖ) given some `labels` a `filename` and refinement factor `ms`."
function create_mesh(
    Lᵢ::Real, Lⱼ::Real, Lₖ::Real,
    labels::Vector,
    filename::String,
    ms::Real=0.5
)

    # Get Labels
    mat_label = labels[1]
    entities_labels = labels[2]
    face_label = entities_labels[1]
    element_label = entities_labels[2]
    bc_labels = labels[3]
    ux_bc_label = bc_labels[1]
    uj_bc_label = bc_labels[2]
    uk_bc_label = bc_labels[3]
    load_label = bc_labels[4]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("$filename")

    # Points
    # Face (x = 0)
    gmsh.model.geo.addPoint(0, 0, 0, ms, 1)
    gmsh.model.geo.addPoint(0, 0, Lₖ, ms, 2)
    gmsh.model.geo.addPoint(0, Lⱼ, Lₖ, ms, 3)
    gmsh.model.geo.addPoint(0, Lⱼ, 0, ms, 4)
    # Face (x = Lᵢ)
    gmsh.model.geo.addPoint(Lᵢ, 0, 0, ms, 5)
    gmsh.model.geo.addPoint(Lᵢ, 0, Lₖ, ms, 6)
    gmsh.model.geo.addPoint(Lᵢ, Lⱼ, Lₖ, ms, 7)
    gmsh.model.geo.addPoint(Lᵢ, Lⱼ, 0, ms, 8)

    # Lines
    gmsh.model.geo.addLine(4, 3, 1)
    gmsh.model.geo.addLine(3, 7, 2)
    gmsh.model.geo.addLine(7, 8, 3)
    gmsh.model.geo.addLine(8, 4, 4)
    gmsh.model.geo.addLine(1, 2, 5)
    gmsh.model.geo.addLine(2, 6, 6)
    gmsh.model.geo.addLine(6, 5, 7)
    gmsh.model.geo.addLine(5, 1, 8)
    gmsh.model.geo.addLine(1, 4, 9)
    gmsh.model.geo.addLine(3, 2, 10)
    gmsh.model.geo.addLine(5, 8, 11)
    gmsh.model.geo.addLine(7, 6, 12)

    # Surfaces
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addCurveLoop([9, -4, -11, 8], 2)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.addCurveLoop([11, -3, 12, 7], 3)
    gmsh.model.geo.addPlaneSurface([3], 3)
    gmsh.model.geo.addCurveLoop([10, -2, -12, 6], 4)
    gmsh.model.geo.addPlaneSurface([4], 4)
    gmsh.model.geo.addCurveLoop([5, -10, -1, -9], 5)
    gmsh.model.geo.addPlaneSurface([5], 5)
    gmsh.model.geo.addCurveLoop([-8, -7, -6, -5], 6)
    gmsh.model.geo.addPlaneSurface([6], 6)
    gmsh.model.geo.addSurfaceLoop([5, 6, 2, 1, 4, 3], 1)

    #Physical names
    dim_surface = 2
    dim_vol = 3

    # Physical Surfaces
    gmsh.model.addPhysicalGroup(dim_surface, [5], 1)
    gmsh.model.addPhysicalGroup(dim_surface, [6], 2)
    gmsh.model.addPhysicalGroup(dim_surface, [2], 3)
    gmsh.model.addPhysicalGroup(dim_surface, [3], 4)

    #Volume
    gmsh.model.geo.addVolume([1], 1)
    gmsh.model.addPhysicalGroup(dim_vol, [1], 5)

    gmsh.model.setPhysicalName(dim_surface, 1, "_$(face_label)_$(ux_bc_label)")
    gmsh.model.setPhysicalName(dim_surface, 2, "_$(face_label)_$(uj_bc_label)")
    gmsh.model.setPhysicalName(dim_surface, 3, "_$(face_label)_$(uk_bc_label)")
    gmsh.model.setPhysicalName(dim_surface, 4, "_$(face_label)_$(load_label)")
    gmsh.model.setPhysicalName(dim_vol, 5, "$(mat_label)_$(element_label)_")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(dim_vol)
    file_name_msh = joinpath(@__DIR__, "./$filename.msh")
    gmsh.write(file_name_msh)
    gmsh.finalize()

    return file_name_msh
end