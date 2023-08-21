"""
Module to export ONSAS data structures to `.vtk` files suitable for [Paraview](https://www.paraview.org/).
"""
module VTK

using Reexport
using WriteVTK

using ..Meshes
using ..Nodes

@reexport import WriteVTK: vtk_grid
export write_vtks

function write_vtks(states_sol, filename)
    nodes_mat = node_matrix(states_sol.analysis.s.mesh)
    connec_mat = connectivity(states_sol.analysis.s.mesh)
    num_elem = length(states_sol.analysis.s.mesh.elements)
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, connec_mat[e])
             for e in 1:num_elem]

    vtk_grid(filename, nodes_mat, cells) do vtk
        vtk["Displacements", VTKPointData()] = rand(3)
    end
end

end
