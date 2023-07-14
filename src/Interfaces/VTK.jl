"""
Module to export ONSAS data structures to `.vtk` files suitable for [Paraview](https://www.paraview.org/).
"""
module VTK

using Reexport
using WriteVTK

using ..Meshes
using ..Nodes
using ..Elements

@reexport import WriteVTK: vtk_grid
export write_vtks

# Add code here.

function node_coordinates_matrix(nodes::Vector{N}) where {dim,T,N<:Node{dim,T}}
    nodes_coords_matrix = Matrix{Float64}(undef, (length(nodes), dim))
    dict_nodes_enum = Dict{AbstractNode,Int}()
    for (i, n) in enumerate(nodes)
        nodes_coords_matrix[i, :] = collect(coordinates(n))
        dict_nodes_enum[n] = i
    end
    return nodes_coords_matrix, dict_nodes_enum
end

function write_vtks(states_sol, filename)
    nodes_mat, nodes_dicts = node_coordinates_matrix(states_sol.analysis.s.mesh.nodes)
    points = transpose(nodes_mat)

    cells = [MeshCell(VTKCellTypes.VTK_TETRA, nodes_dicts.(e.nodes))
             for e in states_sol.analysis.s.mesh.elements]

    vtk_grid(filename, points, cells) do vtk
        # add datasets...
    end
end

end
