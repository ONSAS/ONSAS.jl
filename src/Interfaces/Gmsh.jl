"""
Module defining meshes created in GMSH software.
"""
module Gmsh

using Reexport
using Dictionaries: dictionary
using MshReader: MshFileReader

using ..Entities
using ..Nodes

@reexport import ..Entities: dimension, nodes

# Physical nodes index, all gmsh files should be defined with this label
const PHYSICAL_NODE_LABEL = "node"

export MshFile, connectivity, material_label, entity_label, bc_label, physical_index, gmsh_println,
       PHYSICAL_NODE_LABEL

"""
Collection of `Node`s.
"""
struct MshFile{dim,T,S,I1<:Integer,I2<:Integer}
    "Path to the `.geo` file."
    filename::String
    "Nodes of the mesh."
    vec_nodes::Vector{Node{dim,T}}
    "Connectivity of the mesh."
    connectivity::Vector{Vector{I1}}
    "Physical names."
    physical_index::Vector{I2}
    "Material type labels."
    material_labels::Vector{S}
    "Face and element type labels."
    entities_labels::Vector{S}
    "Boundary condition type labels."
    bc_labels::Vector{S}
    function MshFile(filename::String, vec_nodes::Vector{Node{dim,T}},
                     connectivity::Vector{Vector{I1}}, physical_index::Vector{I2},
                     material_labels::Vector{S}, entities_labels::Vector{S},
                     bc_labels::Vector{S}) where {dim,T,S,I1<:Integer,I2<:Integer}
        @assert length(physical_index) == length(connectivity) "The number of physical indexes = $(length(physical_index))
            must be equal to the number of elements = $(length(connectivity))."
        @assert length(material_labels) == length(entities_labels) == length(bc_labels) "The number of material labels = $(length(material_labels)), 
        entities labels = $(length(entities_labels)) and boundary conditions labels = $(length(bc_labels)) must be equal."

        return new{dim,T,S,I1,I2}(filename, vec_nodes, connectivity, physical_index,
                                  material_labels, entities_labels, bc_labels)
    end
end

"Return material, entities and boundary conditions labels defined in `physical_names`."
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

    return MshFile(filename, vec_nodes, connectivity, physical_index, material_labels,
                   entities_labels, bcs_labels)
end

"Return the connectivity defined in the `MshFile` `mfile`."
connectivity(mfile::MshFile) = mfile.connectivity

"Return physical indexes defined in the `MshFile` `mfile`."
physical_index(mfile::MshFile) = mfile.physical_index

"Return physical_index defined in the `MshFile` `mfile`."
physical_index(mfile::MshFile, entity_index::Integer) = physical_index(mfile)[entity_index]

"Return the mesh dimension defined in the `MshFile` `mfile`."
dimension(::MshFile{dim}) where {dim} = dim

"Return the `Nodes`s of the `MshFile` `mfile`."
nodes(mfile::MshFile) = mfile.vec_nodes

"Return material labels defined in the `MshFile``mfile`."
material_label(mfile::MshFile) = mfile.material_labels

"Return material label for the `AbstractEntity` with index `entity_index`."
function material_label(mfile::MshFile, entity_index::Integer)
    return material_label(mfile)[physical_index(mfile)[entity_index]]
end

"Return entities labels defined in the `MshFile``mfile`."
entity_label(mfile::MshFile) = mfile.entities_labels

"Return entities label for the `AbstractEntity` with index `entity_index`."
function entity_label(mfile::MshFile, entity_index::Integer)
    return entity_label(mfile)[physical_index(mfile)[entity_index]]
end

"Return boundary conditions labels defined in the `MshFile``mfile`."
bc_label(mfile::MshFile) = mfile.bc_labels

"Return boundary conditions label for the `AbstractEntity` with index `entity_index`."
function bc_label(mfile::MshFile, entity_index::Integer)
    return bc_label(mfile)[physical_index(mfile)[entity_index]]
end

"Return the as variables prints from the REPL."
function gmsh_println(output)
    # Extract the number of nodes and elements as integers.
    num_nodes = parse(Int, match(r"(\d+) nodes (\d+) elements", output).captures[1])
    num_entities = parse(Int, match(r"(\d+) nodes (\d+) elements", output).captures[2])
    println("The mesh contains $num_entities entities and $num_nodes nodes. ")
    return dictionary([:nnodes => num_nodes, :num_entities => num_entities])
end

end # module
