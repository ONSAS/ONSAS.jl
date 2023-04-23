using Gmsh

"Creates a mesh for a cylinder loaded with a pressure in its internal face 
and fixed in k at (0,0,0) and (0,0,Lₖ), the boundary conditions and physical 
properties are passed into as `labels`. The .msh file is generated at `filename`
at a given `dir`ectory  with a refinement factor `ms`."
function create_cylinder_mesh(Rᵢ::Real, Rₑ::Real, Lₖ::Real,
                              labels::Vector,
                              filename::String,
                              ms::Real=1,
                              dir=@__DIR__)

    # Refinement factors for internal and external faces
    factorₑ = 0.015
    factorᵢ = factorₑ * Rᵢ / (1.5 * Rₑ)
    msₑ = factorₑ / ms * 2 * pi * Rₑ
    msᵢ = factorᵢ / ms * 2 * pi * Rₑ
    # Get Labels
    mat_label = labels[1]
    entities_labels = labels[2]
    node_label = entities_labels[1]
    face_label = entities_labels[2]
    element_label = entities_labels[3]
    bc_labels = labels[3]
    uᵢ_bc_label = bc_labels[1]
    uⱼ_bc_label = bc_labels[2]
    uₖ_bc_label = bc_labels[3]
    load_label = bc_labels[4]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("$filename")

    # Points

    # Face (z = 0)
    # origin
    gmsh.model.geo.addPoint(0, 0, 0, msᵢ, 1)
    # Rᵢ
    gmsh.model.geo.addPoint.(Rᵢ, 0, 0, msᵢ, 2)
    gmsh.model.geo.addPoint.(0, -Rᵢ, 0, msᵢ, 3)
    gmsh.model.geo.addPoint.(-Rᵢ, 0, 0, msᵢ, 4)
    gmsh.model.geo.addPoint.(0, Rᵢ, 0, msᵢ, 5)
    # Rₑ
    gmsh.model.geo.addPoint.(Rₑ, 0, 0, msₑ, 6)
    gmsh.model.geo.addPoint.(0, -Rₑ, 0, msₑ, 7)
    gmsh.model.geo.addPoint.(-Rₑ, 0, 0, msₑ, 8)
    gmsh.model.geo.addPoint.(0, Rₑ, 0, msₑ, 9)

    # Face (z = Lₖ)
    # origin
    gmsh.model.geo.addPoint.(0, 0, Lₖ, msᵢ, 10)
    # Rᵢ
    gmsh.model.geo.addPoint.(Rᵢ, 0, Lₖ, msᵢ, 11)
    gmsh.model.geo.addPoint.(0, -Rᵢ, Lₖ, msᵢ, 12)
    gmsh.model.geo.addPoint.(-Rᵢ, 0, Lₖ, msᵢ, 13)
    gmsh.model.geo.addPoint.(0, Rᵢ, Lₖ, msᵢ, 14)
    # Rₑ
    gmsh.model.geo.addPoint.(Rₑ, 0, Lₖ, msₑ, 15)
    gmsh.model.geo.addPoint.(0, -Rₑ, Lₖ, msₑ, 16)
    gmsh.model.geo.addPoint.(-Rₑ, 0, Lₖ, msₑ, 17)
    gmsh.model.geo.addPoint.(0, Rₑ, Lₖ, msₑ, 18)

    # Circles (z = 0)
    # Rᵢ
    gmsh.model.geo.addCircleArc(2, 1, 3, 1)
    gmsh.model.geo.addCircleArc(3, 1, 4, 2)
    gmsh.model.geo.addCircleArc(4, 1, 5, 3)
    gmsh.model.geo.addCircleArc(5, 1, 2, 4)
    # Rₑ
    gmsh.model.geo.addCircleArc(6, 1, 7, 5)
    gmsh.model.geo.addCircleArc(7, 1, 8, 6)
    gmsh.model.geo.addCircleArc(8, 1, 9, 7)
    gmsh.model.geo.addCircleArc(9, 1, 6, 8)

    # Circles (z = Lₖ)
    # Rₑ
    gmsh.model.geo.addCircleArc(11, 10, 12, 11)
    gmsh.model.geo.addCircleArc(12, 10, 13, 12)
    gmsh.model.geo.addCircleArc(13, 10, 14, 13)
    gmsh.model.geo.addCircleArc(14, 10, 11, 14)
    # Rᵢ
    gmsh.model.geo.addCircleArc(15, 10, 16, 15)
    gmsh.model.geo.addCircleArc(16, 10, 17, 16)
    gmsh.model.geo.addCircleArc(17, 10, 18, 17)
    gmsh.model.geo.addCircleArc(18, 10, 15, 18)

    # Lines to link face at z = 0 and z = Lₖ
    # Rᵢ
    gmsh.model.geo.addLine(2, 11, 21)
    gmsh.model.geo.addLine(3, 12, 22)
    gmsh.model.geo.addLine(4, 13, 23)
    gmsh.model.geo.addLine(5, 14, 24)
    # Rₑ
    gmsh.model.geo.addLine(6, 15, 25)
    gmsh.model.geo.addLine(7, 16, 26)
    gmsh.model.geo.addLine(8, 17, 27)
    gmsh.model.geo.addLine(9, 18, 28)

    # Curve and surfaces 
    # Rₑ
    gmsh.model.geo.addCurveLoop([-5, 25, 15, -26], 1)
    gmsh.model.geo.addSurfaceFilling([1], 1)
    gmsh.model.geo.addCurveLoop([-6, 26, 16, -27], 2)
    gmsh.model.geo.addSurfaceFilling([2], 2)
    gmsh.model.geo.addCurveLoop([-7, 27, 17, -28], 3)
    gmsh.model.geo.addSurfaceFilling([3], 3)
    gmsh.model.geo.addCurveLoop([-8, 28, 18, -25], 4)
    gmsh.model.geo.addSurfaceFilling([4], 4)
    # Rᵢ
    gmsh.model.geo.addCurveLoop([-1, 21, 11, -22], 5)
    gmsh.model.geo.addSurfaceFilling([-5], 5)
    gmsh.model.geo.addCurveLoop([-2, 22, 12, -23], 6)
    gmsh.model.geo.addSurfaceFilling([-6], 6)
    gmsh.model.geo.addCurveLoop([-3, 23, 13, -24], 7)
    gmsh.model.geo.addSurfaceFilling([-7], 7)
    gmsh.model.geo.addCurveLoop([-4, 24, 14, -21], 8)
    gmsh.model.geo.addSurfaceFilling([-8], 8)

    # Surface at z = 0 
    gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 9)
    gmsh.model.geo.addCurveLoop([-4, -3, -2, -1], 10)
    gmsh.model.geo.addPlaneSurface([9, 10], 9)

    # Surface at z = Lₖ 
    gmsh.model.geo.addCurveLoop([-18, -17, -16, -15], 11)
    gmsh.model.geo.addCurveLoop([11, 12, 13, 14], 12)
    gmsh.model.geo.addPlaneSurface([11, 12], 10)

    # Volume surface
    gmsh.model.geo.addSurfaceLoop([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 1)
    gmsh.model.geo.addVolume([1], 1)

    # Set physical to apply boundary conditions
    #Physical names
    dim_node = 0
    dim_surface = 2
    dim_vol = 3

    # nodes
    # fixed uᵢ
    gmsh.model.addPhysicalGroup(dim_node, [9, 18], 1)
    gmsh.model.setPhysicalName(dim_node, 1, "_$(node_label)_$(uᵢ_bc_label)")
    # fixed uⱼ
    gmsh.model.addPhysicalGroup(dim_node, [6, 15], 2)
    gmsh.model.setPhysicalName(dim_node, 2, "_$(node_label)_$(uⱼ_bc_label)")
    # surfaces
    # fixed uₖ
    gmsh.model.addPhysicalGroup(dim_surface, [9, 10], 3)
    gmsh.model.setPhysicalName(dim_surface, 3, "_$(face_label)_$(uₖ_bc_label)")

    # pressure
    gmsh.model.addPhysicalGroup(dim_surface, [5, 6, 7, 8], 4)
    gmsh.model.setPhysicalName(dim_surface, 4, "_$(face_label)_$(load_label)")

    # volume
    gmsh.model.addPhysicalGroup(dim_vol, [1], 5)
    gmsh.model.setPhysicalName(dim_vol, 5, "$(mat_label)_$(element_label)_")

    # synchronize and write 
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(dim_vol)
    filename_msh = joinpath(dir, filename * ".msh")
    gmsh.write(filename_msh)
    gmsh.finalize()

    return filename_msh
end
