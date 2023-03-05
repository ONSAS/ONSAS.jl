"""
Module defining meshes created in GMSH software.
"""
module Gmsh

using MshReader: MshFileReader
using ..Elements: Node

export MSHFile


""" MSHFile.
A `MSHFile` is a collection of `Node`s. Also the filename   
### Fields:
- `filename`              -- Stores the .geo file name.
- `vec_nodes`              -- Stores the `Nodes`s of the mesh.
- `connectivity`          -- Stores a `Vector` of `Vector` with the node indexes.
- `material_labels`       -- Stores a `Vector` of `String`s with the material type labels defined in the .geo.
- `element_face_labels`   -- Stores a `Vector` of `String`s with the face and element type labels defined in the .geo.
- `bc_labels`             -- Stores a `Vector` of `String`s with the boundary condition type labels defined in the .geo.
"""
struct MSHFile{dim,T,S,I1<:Integer,I2<:Integer}
    filename::String
    vec_nodes::Vector{Node{dim,T}}
    connectivity::Vector{Vector{I1}}
    physical_indexes::Vector{I2}
    material_labels::Vector{S}
    entities_labels::Vector{S}
    bc_labels::Vector{S}
end

"Constructor of the `MSHFile` object with a file name `filename`."
function MSHFile(filename::String)
    nodes_coords, connectivity, physical_names, physical_indexes = MshFileReader(filename)

    material_labels, entities_labels, bcs_labels = _getlabels(physical_names)

    vec_nodes = [Node(n) for n in eachrow(nodes_coords)]

    MSHFile(filename, vec_nodes, connectivity, physical_indexes, material_labels, entities_labels, bcs_labels)
end

function _getlabels(physical_names::Vector{String})

    physical_names = split.(physical_names, "_")

    # Entities material, types and boundary conditions
    length_labels_elements = 3
    material_labels = getindex.(filter(l -> length(l) == length_labels_elements, physical_names), 1)
    entity_labels = getindex.(filter(l -> length(l) == length_labels_elements, physical_names), 2)
    bcs_labels = getindex.(filter(l -> length(l) == length_labels_elements, physical_names), 3)

    return material_labels, entity_labels, bcs_labels

end

end #endModule