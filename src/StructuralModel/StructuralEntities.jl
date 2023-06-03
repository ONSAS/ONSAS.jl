using Dictionaries
using Reexport

using ..Elements
using ..Utils
using ..Meshes
using ..Gmsh

@reexport import ..Meshes: Mesh

export StructuralEntities, all_entities, face_types_to_faces, face_types, elem_types_to_elements,
       elem_types

"""
Struct used to define meshes via GMSH.

TODO Improve description. Parameter for inner type of VE?
Consists of a collection of element and face types, each assigned to a corrsponding vector of `Face`s and `Element`s. 
"""
struct StructuralEntities{F<:AbstractFace,VF<:Vector,E<:AbstractElement,VE<:Vector}
    "Dictionary with `Element` types (`Element`s without assigned `Node`s) as keys and the corresponding elements created."
    elem_types_to_elements::Dictionary{E,VE}
    "Dictionary with `Face` types (`Face`s without assigned `Node`s) as keys and the corresponding faces created."
    face_types_to_faces::Dictionary{F,VF}
    function StructuralEntities(elem_types_to_elements::Dictionary{E,VE},
                                face_types_to_faces::Dictionary{F,VF}) where
             {F<:AbstractFace,VF<:Vector,E<:AbstractElement,VE<:Vector}
        velems = collect(keys(elem_types_to_elements))
        vfaces = collect(keys(face_types_to_faces))
        vlabels = vcat(label.(velems), label.(vfaces)) # mapreduce
        @assert length(vlabels) == length(unique(vlabels)) "Every `Face` and `Element` type labels must be different" # allunique
        return new{F,VF,E,VE}(elem_types_to_elements, face_types_to_faces)
    end
end

"Constructor for an empty `StructuralEntities` with a `Vector` of `Element`s `velems` and `Face`s `vfaces`."
function StructuralEntities(velems::Vector{E},
                            vfaces::Vector{F}=Vector{AbstractFace}()) where {E<:AbstractElement,
                                                                             F<:AbstractFace}
    elem_types_to_elements = dictionary(map(elem -> elem => Vector{typeof(elem)}(), velems))
    face_types_to_faces = dictionary(map(face -> face => Vector{typeof(face)}(), vfaces))
    return StructuralEntities(elem_types_to_elements, face_types_to_faces)
end

"Return a `Dictionary` with `Element` types as keys and the corresponding `Element`s as values."
elem_types_to_elements(entities::StructuralEntities) = entities.elem_types_to_elements

"Return a `Dictionary` with `Face` types as keys and the corresponding `Face`s as values."
face_types_to_faces(entities::StructuralEntities) = entities.face_types_to_faces

"Return the `Vector` of `Element` types defined in the `StructuralEntities` `entities`."
elem_types(entities::StructuralEntities) = collect(keys(entities.elem_types_to_elements))

"Return the `Vector` of `Face` types defined in the `StructuralEntities` `entities`."
face_types(entities::StructuralEntities) = collect(keys(entities.face_types_to_faces))

"Return all `Entity`s defined into `StructuralEntities`."
function all_entities(entities::StructuralEntities)
    return unique(vcat(face_types(entities), elem_types(entities))) # unique! ?
end

"Return the `Entity` with the label `l` in the `StructuralEntities` `entities`."
function Base.getindex(entities::StructuralEntities, l::L) where {L<:Union{Symbol,AbstractString}}
    return first(filter(ent -> label(ent) == Symbol(l), all_entities(entities)))
end

"Return the `Mesh` given an `MshFile` `msh_file` and a `StructuralEntities` `entities`."
function Mesh(msh_file::MshFile, entities::StructuralEntities)

    # Initialize the mesh with its nodes
    nodes = Meshes.nodes(msh_file)
    mesh = Mesh(; nodes, extra=msh_file)

    # Loop over all physical entities and push them into the mesh
    for (entity_index, entity_nodes_indexes) in enumerate(connectivity(msh_file))

        # Create entity and push it into the mesh
        nodes_entity = view(nodes, entity_nodes_indexes)
        entity_type_label = entity_label(msh_file, entity_index)

        # Check if the entity is a node, if not add it to the mesh
        local entity_position # in the mesh vetor of of entities

        if entity_type_label == PHYSICAL_NODE_LABEL
            entity_position = entity_nodes_indexes[]
            entity = nodes_entity[]
        else
            entity_type = entities[entity_type_label]
            entity = create_entity(entity_type, nodes_entity)
            push!(mesh, entity)
            entity_position = length(mesh, entity)
        end
        add_entity_to_set!(mesh, String(entity_type_label), entity_position, entity)

        # Find material and push into the set
        material_type_label = material_label(msh_file, entity_index)
        if !isempty(material_type_label)
            add_entity_to_set!(mesh, String(material_type_label), entity_position, entity)
        end

        # Find boundary conditions 
        bc_type_label = bc_label(msh_file, entity_index)
        if !isempty(bc_type_label)
            add_entity_to_set!(mesh, String(bc_type_label), entity_position, entity)
        end
    end
    return mesh
end
