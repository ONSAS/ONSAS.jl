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

    n_times = length(displacements(states_sol))
    n_dofs = length(displacements(states_sol)[end])
    mypad = Integer(ceil(log10(n_times))) + 1
    for i in 1:n_times
        pdata = displacements(states_sol)[i]
        nodes = nodes_mat + reshape(pdata, (3, Integer(n_dofs / 3))) # TO DO generalize for angle dof cases
        filename_i = filename * string(i; pad=mypad)
        vtk_grid(filename_i, nodes, cells) do vtk
            vtk["Displacements", VTKPointData()] = pdata
        end
    end
end

end
