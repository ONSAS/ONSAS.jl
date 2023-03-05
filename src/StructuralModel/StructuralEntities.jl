using ..Elements: AbstractElement, AbstractFace
using Dictionaries: Dictionary, dictionary
using ..Utils: label
using ..Meshes: MSHFile

import ..Meshes: Mesh

export StructuralEntities

""" Structural elements struct.
A `StructuralMaterials` is a collection of `Element`s and `Faces`s types assigning to a vector of `Face`s and `Element`s.
This struct is used to define meshes via GMSH. 
### Fields:
- `elem_types_to_elements` -- Store a dictionary with `Element` types (`Element`s without assigned `Node`s) 
as keys and the corresponding elements created.
- `face_types_to_faces` -- Store a dictionary with `Face` types (`Face`s without assigned `Node`s) 
as keys and the corresponding faces created.
"""

struct StructuralEntities{F<:AbstractFace,E<:AbstractElement}
    elem_types_to_elements::Dictionary{E,Vector{E}}
    face_types_to_faces::Dictionary{F,Vector{F}}
    function StructuralEntities(elem_types_to_elements::Dictionary{E,Vector{E}}, face_types_to_faces::Dictionary{F,Vector{F}}) where {F<:AbstractFace,E<:AbstractElement}
        velems = collect(keys(elem_types_to_elements))
        vfaces = collect(keys(face_types_to_faces))
        vlabels = vcat(label.(velems), label.(vfaces))
        @assert length(vlabels) == length(unique(vlabels)) "Every `Face` and `Element` type labels must be different"
        new{F,E}(elem_types_to_elements, face_types_to_faces)
    end
end


"Constructor for an empty `StructuralEntities` with a `Vector` of `Element`s `velems` and `Face`s `vfaces`."
function StructuralEntities(velems::Vector{E}, vfaces::Vector{F}=Vector{AbstractFace}()) where {E<:AbstractElement,F<:AbstractFace}
    elem_types_to_elements = dictionary(map(elem -> elem => Vector{typeof(elem)}(), velems))
    face_types_to_faces = dictionary(map(face -> face => Vector{typeof(face)}(), vfaces))
    StructuralEntities(elem_types_to_elements, face_types_to_faces)
end


function Mesh(msh_file::MSHFile, s_entities::StructuralEntities)



    for (index_entity, entity_nodes_indexes) in enumerate(eachcol(msh_file.connectivity))
        nodes_entity = [msh_file.vec_nodes[node_index] for node_index in entity_nodes_indexes]
        entity_type_label = mesh_file.entities_labels[index_entity]
        entity_type = s_entities[entity_type_label]
        entity = create_entity(entity_type, nodes_entity)

end