using Gmsh

"Creates a mesh for a cube with a loaded face and a fixed origin."
function create_linear_extension_mesh(Lᵢ::Real, Lⱼ::Real, Lₖ::Real,
                                      labels::Vector,
                                      filename::String,
                                      ms::Real=0.5,
                                      dir::String=@__DIR__)

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
    #num node1, num node 2, num line 
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
    # Top surface
    top_surface_index = 1
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], top_surface_index)
    gmsh.model.geo.addPlaneSurface([1], top_surface_index)
    # Front surface
    front_surface_index = 2
    gmsh.model.geo.addCurveLoop([9, -4, -11, 8], front_surface_index)
    gmsh.model.geo.addPlaneSurface([2], front_surface_index)
    # Right surface
    right_surface_index = 3
    gmsh.model.geo.addCurveLoop([11, -3, 12, 7], right_surface_index)
    gmsh.model.geo.addPlaneSurface([3], right_surface_index)
    # Back surface
    back_surface_index = 4
    gmsh.model.geo.addCurveLoop([10, -2, -12, 6], back_surface_index)
    gmsh.model.geo.addPlaneSurface([4], back_surface_index)
    # Left surface
    left_surface_index = 5
    gmsh.model.geo.addCurveLoop([5, -10, -1, -9], left_surface_index)
    gmsh.model.geo.addPlaneSurface([5], left_surface_index)
    # Right surface
    bottom_surface_index = 6
    gmsh.model.geo.addCurveLoop([-8, -7, -6, -5], bottom_surface_index)
    gmsh.model.geo.addPlaneSurface([6], bottom_surface_index)
    # Volume
    gmsh.model.geo.addSurfaceLoop([5, 6, 2, 1, 4, 3], 1)

    #Physical names
    dim_surface = 2
    dim_vol = 3

    # Physical Surfaces
    fixed_ux_surfaces_index = 1
    gmsh.model.addPhysicalGroup(dim_surface, [left_surface_index], fixed_ux_surfaces_index)
    fixed_uy_surfaces_index = 2
    gmsh.model.addPhysicalGroup(dim_surface, [bottom_surface_index, top_surface_index],
                                fixed_uy_surfaces_index)
    fixed_uz_surfaces_index = 3
    gmsh.model.addPhysicalGroup(dim_surface, [front_surface_index, back_surface_index],
                                fixed_uz_surfaces_index)
    load_surfaces_index = 4
    gmsh.model.addPhysicalGroup(dim_surface, [right_surface_index], load_surfaces_index)

    #Volume
    gmsh.model.geo.addVolume([1], 1)
    gmsh.model.addPhysicalGroup(dim_vol, [1], 5)

    gmsh.model.setPhysicalName(dim_surface, fixed_ux_surfaces_index,
                               "_$(face_label)_$(ux_bc_label)")
    gmsh.model.setPhysicalName(dim_surface, fixed_uy_surfaces_index,
                               "_$(face_label)_$(uj_bc_label)")
    gmsh.model.setPhysicalName(dim_surface, fixed_uz_surfaces_index,
                               "_$(face_label)_$(uk_bc_label)")
    gmsh.model.setPhysicalName(dim_surface, load_surfaces_index, "_$(face_label)_$(load_label)")
    gmsh.model.setPhysicalName(dim_vol, 5, "$(mat_label)_$(element_label)_")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(dim_vol)
    filename_msh = joinpath(dir, filename * ".msh")
    gmsh.write(filename_msh)
    gmsh.finalize()

    filename_msh
end
