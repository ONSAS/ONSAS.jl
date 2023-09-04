"""
Module to export ONSAS data structures to `.vtk` files suitable for [Paraview](https://www.paraview.org/).
"""
module VTK

using Reexport
using WriteVTK

using ..Meshes
using ..Nodes
using ..Solutions
using ..Structures
using ..StructuralAnalyses

export write_vtk
@reexport import WriteVTK: vtk_grid, write_vtk

"""
Generate a VTK file given a solution struct.
Currently only `displacements` are exported for each time point.
"""
function write_vtk(sol::AbstractSolution, filename::String)
    msh = mesh(structure(analysis(sol)))
    nodes_mat = node_matrix(msh)
    connec_mat = connectivity(msh)
    num_elem = length(elements(msh))
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, connec_mat[e]) for e in 1:num_elem]

    n_times = length(displacements(sol))
    n_dofs = num_dofs(msh)
    mypad = Integer(ceil(log10(n_times))) + 1
    for i in 1:n_times
        pdata = displacements(sol)[i]
        nodes = nodes_mat + reshape(pdata, (3, Integer(n_dofs / 3))) # TO DO generalize for angle dof cases
        filename_i = filename * string(i; pad=mypad)
        vtk_grid(filename_i, nodes, cells) do vtk
            vtk["Displacements", VTKPointData()] = pdata
        end
    end
    @debug "VTK output written to $filename"
    filename
end

end
