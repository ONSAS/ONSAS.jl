"""
Module to export ONSAS data structures to `.vtk` files suitable for [Paraview](https://www.paraview.org/).
"""
module VTK

using Reexport
using WriteVTK

using ..Meshes
using ..Nodes
using ..StructuralAnalyses

@reexport import WriteVTK: vtk_grid
export write_vtks

function write_vtks(states_sol, filename)
    nodes_mat = node_matrix(states_sol.analysis.s.mesh)
    connec_mat = connectivity(states_sol.analysis.s.mesh)
    num_elem = length(states_sol.analysis.s.mesh.elements)
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, connec_mat[e]) for e in 1:num_elem]

    s = states_sol.analysis.s
    nnodes = length(nodes(s))
    vtk_grid(filename, nodes_mat, cells) do vtk
        vtk["Displacements x", VTKPointData()] = rand(3 * nnodes)
    end

    # Main.@infiltrate
    #=
    s = states_sol.analysis.s
    for n in nodes(s)
        # Get the displacements for all coordinates at the last time point.
        data = [last(c) for c in displacements(states_sol, n)]
        vtk_point_data(filename, data, "Displacements")
    end
    =#
end

end
