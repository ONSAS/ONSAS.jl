"""
Module defining meshes created in GMSH software.
"""
module Gmsh

using Reexport: @reexport
using MshReader: MshFileReader
using ..Elements: Node

@reexport import ..Elements: dimension, nodes

export MshFile, connectivity, material_label, entity_label, bc_label, physical_index


""" MshFile.
A `MshFile` is a collection of `Node`s. Also the filename   
### Fields:
- `filename`              -- Stores the .geo file name.
- `vec_nodes`              -- Stores the `Nodes`s of the mesh.
- `connectivity`          -- Stores a `Vector` of `Vector` with the node indexes.
- `material_labels`       -- Stores a `Vector` of `String`s with the material type labels defined in the .geo.
- `element_face_labels`   -- Stores a `Vector` of `String`s with the face and element type labels defined in the .geo.
- `bc_labels`             -- Stores a `Vector` of `String`s with the boundary condition type labels defined in the .geo.
"""
struct MshFile{dim,T,S,I1<:Integer,I2<:Integer}
    filename::String
    vec_nodes::Vector{Node{dim,T}}
    connectivity::Vector{Vector{I1}}
    physical_index::Vector{I2}
    material_labels::Vector{S}
    entities_labels::Vector{S}
    bc_labels::Vector{S}
    function MshFile(filename::String, vec_nodes::Vector{Node{dim,T}},
        connectivity::Vector{Vector{I1}}, physical_index::Vector{I2},
        material_labels::Vector{S}, entities_labels::Vector{S}, bc_labels::Vector{S}
    ) where {dim,T,S,I1<:Integer,I2<:Integer}

        @assert length(physical_index) == length(connectivity) "The number of physical indexes = $(length(physical_index))
            must be equal to the number of elements = $(length(connectivity))."
        @assert length(material_labels) == length(entities_labels) == length(bc_labels) "The number of material labels = $(length(material_labels)), 
        entities labels = $(length(entities_labels)) and boundary conditions labels = $(length(bc_labels)) must be equal."

        return new{dim,T,S,I1,I2}(filename, vec_nodes, connectivity, physical_index, material_labels, entities_labels, bc_labels)
    end
end

"Returns material, entities and boundary conditions labels defined in `physical_names`."
function _getlabels(physical_names::Vector{String})

    split_char = "_"

    physical_names = split.(physical_names, split_char)

    # Entities material, types and boundary conditions
    length_labels_elements = 3
    material_labels = getindex.(filter(l -> length(l) == length_labels_elements, physical_names), 1)
    entity_labels = getindex.(filter(l -> length(l) == length_labels_elements, physical_names), 2)
    bcs_labels = getindex.(filter(l -> length(l) == length_labels_elements, physical_names), 3)

    return material_labels, entity_labels, bcs_labels

end

"Constructor of the `MshFile` object with a file name `filename`."
function MshFile(filename::String)
    nodes_coords, connectivity, physical_names, physical_index = MshFileReader(filename)

    material_labels, entities_labels, bcs_labels = _getlabels(physical_names)

    vec_nodes = [Node(n) for n in eachrow(nodes_coords)]

    MshFile(filename, vec_nodes, connectivity, physical_index, material_labels, entities_labels, bcs_labels)
end

"Returns the connectivity defined in the `MshFile` `mfile`."
connectivity(mfile::MshFile) = mfile.connectivity

"Returns physical indexes defined in the `MshFile` `mfile`."
physical_index(mfile::MshFile) = mfile.physical_index

"Returns physical_index defined in the `MshFile` `mfile`."
physical_index(mfile::MshFile, entity_index::Integer) = physical_index(mfile)[entity_index]

"Returns the mesh dimension defined in the `MshFile` `mfile`."
dimension(::MshFile{dim}) where {dim} = dim

"Returns the `Nodes`s of the `MshFile` `mfile`."
nodes(mfile::MshFile) = mfile.vec_nodes

"Returns material labels defined in the `MshFile``mfile`."
material_label(mfile::MshFile) = mfile.material_labels

"Returns material label for the `AbstractEntity` with index `entity_index`."
material_label(mfile::MshFile, entity_index::Integer) = material_label(mfile)[physical_index(mfile)[entity_index]]

"Returns entities labels defined in the `MshFile``mfile`."
entity_label(mfile::MshFile) = mfile.entities_labels

"Returns entities label for the `AbstractEntity` with index `entity_index`."
entity_label(mfile::MshFile, entity_index::Integer) = entity_label(mfile)[physical_index(mfile)[entity_index]]

"Returns boundary conditions labels defined in the `MshFile``mfile`."
bc_label(mfile::MshFile) = mfile.bc_labels

"Returns boundary conditions label for the `AbstractEntity` with index `entity_index`."
bc_label(mfile::MshFile, entity_index::Integer) = bc_label(mfile)[physical_index(mfile)[entity_index]]

end # module