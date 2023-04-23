using Dictionaries
using Reexport

using ..Elements
using ..Utils
using ..Meshes
using ..Meshes.Gmsh: PHYSICAL_NODE_LABEL

@reexport import ..Meshes: Mesh

export StructuralEntities, all_entities, face_types_to_faces, face_types, elem_types_to_elements, elem_types

""" Structural elements struct.
A `StructuralMaterials` is a collection of `Element`s and `Faces`s types assigning to a vector of `Face`s and `Element`s.
This struct is used to define meshes via GMSH. 
### Fields:
- `elem_types_to_elements` -- Store a dictionary with `Element` types (`Element`s without assigned `Node`s) 
as keys and the corresponding elements created.
- `face_types_to_faces` -- Store a dictionary with `Face` types (`Face`s without assigned `Node`s) 
as keys and the corresponding faces created.
"""
struct StructuralEntities{F<:AbstractFace,VF<:Vector,E<:AbstractElement,VE<:Vector}
    elem_types_to_elements::Dictionary{E,VE}
    face_types_to_faces::Dictionary{F,VF}
    function StructuralEntities(elem_types_to_elements::Dictionary{E,VE}, face_types_to_faces::Dictionary{F,VF}) where
    {F<:AbstractFace,VF<:Vector,E<:AbstractElement,VE<:Vector}
        velems = collect(keys(elem_types_to_elements))
        vfaces = collect(keys(face_types_to_faces))
        vlabels = vcat(label.(velems), label.(vfaces))
        @assert length(vlabels) == length(unique(vlabels)) "Every `Face` and `Element` type labels must be different"
        new{F,VF,E,VE}(elem_types_to_elements, face_types_to_faces)
    end
end

"Constructor for an empty `StructuralEntities` with a `Vector` of `Element`s `velems` and `Face`s `vfaces`."
function StructuralEntities(velems::Vector{E}, vfaces::Vector{F}=Vector{AbstractFace}()) where {E<:AbstractElement,F<:AbstractFace}
    elem_types_to_elements = dictionary(map(elem -> elem => Vector{typeof(elem)}(), velems))
    face_types_to_faces = dictionary(map(face -> face => Vector{typeof(face)}(), vfaces))
    StructuralEntities(elem_types_to_elements, face_types_to_faces)
end

"Return a `Dictionary` with `Element` types as keys and the corresponding `Element`s as values."
elem_types_to_elements(s_entities::StructuralEntities) = s_entities.elem_types_to_elements

"Return a `Dictionary` with `Face` types as keys and the corresponding `Face`s as values."
face_types_to_faces(s_entities::StructuralEntities) = s_entities.face_types_to_faces

"Return the `Vector` of `Element` types defined in the `StructuralEntities` `s_entities`."
elem_types(s_entities::StructuralEntities) = collect(keys(s_entities.elem_types_to_elements))

"Return the `Vector` of `Face` types defined in the `StructuralEntities` `s_entities`."
face_types(s_entities::StructuralEntities) = collect(keys(s_entities.face_types_to_faces))

"Return all `Entity`s defined into `StructuralEntities`."
all_entities(s_entities::StructuralEntities) = unique(vcat(face_types(s_entities), elem_types(s_entities)))

"Return the `Entity` with the label `l` in the `StructuralEntities` `s_entities`."
Base.getindex(s_entities::StructuralEntities, l::L) where {L<:Union{Symbol,AbstractString}} =
    first(filter(ent -> label(ent) == Symbol(l), all_entities(s_entities)))

"Return the `Mesh` given an `MshFile` `msh_file` and a `StructuralEntities` `s_entities`."
function Mesh(msh_file::MshFile, s_entities::StructuralEntities)

    # Initialize the mesh with its nodes
    nodes = Meshes.nodes(msh_file)
    mesh = Mesh(nodes, msh_file)

    # Loop over all physical entities and push them into the mesh
    for (entity_index, entity_nodes_indexes) in enumerate(connectivity(msh_file))

        # Create entity and push it into the mesh
        nodes_entity = view(nodes, entity_nodes_indexes)
        entity_type_label = entity_label(msh_file, entity_index)

        # Check if the entity is a node, if not add it to the mesh
        local entity_position # in the mesh vetor of of entities

        if entity_type_label == PHYSICAL_NODE_LABEL
            # Since is a node 
            entity_position = entity_nodes_indexes[]
            entity = nodes_entity[]
            _add_entity_to_set!(mesh, entity_type_label, entity_position, entity)
        else
            entity_type = s_entities[entity_type_label]
            entity = create_entity(entity_type, nodes_entity)
            entity_position = push!(mesh, entity)
            # Add entity label into sets
            _add_entity_to_set!(mesh, entity_type_label, entity_position, entity)
        end

        # Find material and push into the set
        material_type_label = material_label(msh_file, entity_index)
        ~isempty(material_type_label) &&
            _add_entity_to_set!(mesh, material_type_label, entity_position, entity)

        # Find boundary conditions 
        bc_type_label = bc_label(msh_file, entity_index)
        ~isempty(bc_type_label) &&
            _add_entity_to_set!(mesh, bc_type_label, entity_position, entity)
    end
    mesh
end

"Add an entity index of type `AbstractNode` to the `Mesh` `m` node `Set`s."
function _add_entity_to_set!(m::Mesh, entity_type_label::S, entity_position::Int, ::AbstractNode) where {S}
    add_node_to_set!(m, Symbol(entity_type_label), entity_position)
end

"Add an entity index of type `AbstractFace` to the `Mesh` `m` face `Set`s."
_add_entity_to_set!(m::Mesh, entity_type_label::S, entity_position::Int, ::AbstractFace) where {S} =
    add_face_to_set!(m, Symbol(entity_type_label), entity_position)

"Add an entity index of type `AbstractElement` to the `Mesh` `m` element `Set`s."
_add_entity_to_set!(mesh::Mesh, entity_type_label::S, entity_position::Int, ::AbstractElement) where {S} =
    add_element_to_set!(mesh, Symbol(entity_type_label), entity_position)
