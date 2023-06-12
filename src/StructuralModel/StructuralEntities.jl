"""
Module defining structural entities properties to then construct an structure. This used
to create a mesh with GMSH, and the structural entities contain the ONSAS element type to
the posterior conversion to the ONSAS mesh type.
"""
module StructuralEntities

using Dictionaries
using Reexport

using ..Utils
using ..Nodes
using ..Entities
using ..Meshes
using ..Gmsh

@reexport import ..Meshes: Mesh

export StructuralEntity, all_entities, face_types_to_faces, face_types, elem_types_to_elements,
       elem_types

"""
Struct used to define meshes via GMSH.
Consists of a collection of element and face types, each assigned to a corresponding vector of `Face`s and `Element`s.
"""
struct StructuralEntity{F<:AbstractFace,VF<:Vector,E<:AbstractElement,VE<:Vector}
    "Dictionary with `Element` types (`Element`s without assigned `Node`s) as keys and the corresponding elements created."
    elem_types_to_elements::Dictionary{E,VE}
    "Dictionary with `Face` types (`Face`s without assigned `Node`s) as keys and the corresponding faces created."
    face_types_to_faces::Dictionary{F,VF}
    function StructuralEntity(elem_types_to_elements::Dictionary{E,VE},
                              face_types_to_faces::Dictionary{F,VF}) where
             {F<:AbstractFace,VF<:Vector,E<:AbstractElement,VE<:Vector}
        velems = collect(keys(elem_types_to_elements))
        vfaces = collect(keys(face_types_to_faces))
        vlabels = vcat(label.(velems), label.(vfaces)) # mapreduce
        @assert length(vlabels) == length(unique(vlabels)) "Every `Face` and `Element` type labels must be different" # allunique
        new{F,VF,E,VE}(elem_types_to_elements, face_types_to_faces)
    end
end

"Constructor for an empty `StructuralEntity` with a `Vector` of `Element`s `velems` and `Face`s `vfaces`."
function StructuralEntity(velems::Vector{E},
                          vfaces::Vector{F}=Vector{AbstractFace}()) where {E<:AbstractElement,
                                                                           F<:AbstractFace}
    elem_types_to_elements = dictionary(map(elem -> elem => Vector{typeof(elem)}(), velems))
    face_types_to_faces = dictionary(map(face -> face => Vector{typeof(face)}(), vfaces))
    StructuralEntity(elem_types_to_elements, face_types_to_faces)
end

"Return a `Dictionary` with `Element` types as keys and the corresponding `Element`s as values."
elem_types_to_elements(entities::StructuralEntity) = entities.elem_types_to_elements

"Return a `Dictionary` with `Face` types as keys and the corresponding `Face`s as values."
face_types_to_faces(entities::StructuralEntity) = entities.face_types_to_faces

"Return the `Vector` of `Element` types defined in the `StructuralEntity` `entities`."
elem_types(entities::StructuralEntity) = collect(keys(entities.elem_types_to_elements))

"Return the `Vector` of `Face` types defined in the `StructuralEntity` `entities`."
face_types(entities::StructuralEntity) = collect(keys(entities.face_types_to_faces))

"Return all `Entity`s defined into `StructuralEntity`."
function all_entities(entities::StructuralEntity)
    unique(vcat(face_types(entities), elem_types(entities))) # unique! ?
end

"Return the `Entity` with the label `l` in the `StructuralEntity` `entities`."
function Base.getindex(entities::StructuralEntity, l::L) where {L<:Union{Symbol,AbstractString}}
    first(filter(ent -> label(ent) == Symbol(l), all_entities(entities)))
end

"Return the `Mesh` given an `MshFile` `msh_file` and a `StructuralEntity` `entities`."
function Mesh(msh_file::MshFile, entities::StructuralEntity)

    # Initialize the mesh vectors
    nodes = Meshes.nodes(msh_file)
    faces = AbstractFace[]
    elements = AbstractElement[]

    # Initialize the mesh sets
    node_set = EntitySet()
    face_set = EntitySet()
    element_set = EntitySet()

    # Loop over all physical entities and push them into the mesh
    for (entity_index, entity_nodes_indexes) in enumerate(connectivity(msh_file))

        # Create entity and push it into the mesh
        nodes_entity = view(nodes, entity_nodes_indexes)
        entity_type_label = entity_label(msh_file, entity_index)

        # Check if the entity is a node, if not add it to the mesh
        local entity_position # in the mesh vetor of of entities
        local entity # type of the entity

        if entity_type_label == PHYSICAL_NODE_LABEL
            entity_position = entity_nodes_indexes[]
            entity = nodes_entity[]
            entity_position
        else
            entity_type = entities[entity_type_label]
            entity = create_entity(entity_type, nodes_entity)
            if entity_type isa AbstractFace
                push!(faces, entity)
                entity_position = length(faces)
            elseif entity_type isa AbstractElement
                push!(elements, entity)
                entity_position = length(elements)
            else
                throw(ArgumentError("Entity type $entity_type not recognized"))
            end
        end

        if entity isa AbstractNode
            # Add entity to set
            add_node_to_set!(node_set, entity_type_label, entity_position)

            # Add boundary condition to set
            bc_type_label = bc_label(msh_file, entity_index)
            add_node_to_set!(node_set, bc_type_label, entity_position)

        elseif entity isa AbstractFace

            # Add entity to set
            add_face_to_set!(face_set, entity_type_label, entity_position)

            # Add boundary condition to set
            bc_type_label = bc_label(msh_file, entity_index)
            add_face_to_set!(face_set, bc_type_label, entity_position)

        elseif entity isa AbstractElement

            # Add entity to set
            add_element_to_set!(element_set, entity_type_label, entity_position)

            # Add boundary condition to set
            bc_type_label = bc_label(msh_file, entity_index)
            add_element_to_set!(element_set, bc_type_label, entity_position)

            # Add material label to set
            material_type_label = material_label(msh_file, entity_index)
            add_element_to_set!(element_set, material_type_label, entity_position)
        end
    end

    # Remove "" labels
    empty_label = ""
    haskey(element_set, empty_label) && delete!(element_set, empty_label)
    haskey(face_set, empty_label) && delete!(face_set, empty_label)
    haskey(node_set, empty_label) && delete!(node_set, empty_label)

    Mesh(nodes, elements, faces, node_set, element_set, face_set, msh_file)
end

end # module
