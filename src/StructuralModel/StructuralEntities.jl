using ..Elements: AbstractElement, AbstractFace
using Dictionaries: Dictionary, dictionary
using ..Utils: label
using ..Meshes: MshFile

import ..Meshes: Mesh

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

"Returns a `Dictionary` with `Element` types as keys and the corresponding `Element`s as values."
elem_types_to_elements(s_entities::StructuralEntities) = s_entities.elem_types_to_elements

"Returns a `Dictionary` with `Face` types as keys and the corresponding `Face`s as values."
face_types_to_faces(s_entities::StructuralEntities) = s_entities.face_types_to_faces

"Returns the `Vector` of `Element` types defined in the `StructuralEntities` `s_entities`."
elem_types(s_entities::StructuralEntities) = collect(keys(s_entities.elem_types_to_elements))

"Returns the `Vector` of `Face` types defined in the `StructuralEntities` `s_entities`."
face_types(s_entities::StructuralEntities) = collect(keys(s_entities.face_types_to_faces))

"Returns all `Entity`s defined into `StructuralEntities`."
all_entities(s_entities::StructuralEntities) = unique(vcat(face_types(s_entities), elem_types(s_entities)))

"Returns the `Entity` with the label `l` in the `StructuralEntities` `s_entities`."
Base.getindex(s_entities::StructuralEntities, l::L) where {L<:Union{Symbol,AbstractString}} =
    first(filter(ent -> label(ent) == Symbol(l), all_entities(s_entities)))
