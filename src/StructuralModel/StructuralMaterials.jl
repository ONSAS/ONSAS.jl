using ..Materials: AbstractMaterial
using ..Elements: AbstractElement
using Dictionaries: Dictionary, dictionary
using ..Utils: label

export StructuralMaterials

""" Structural materials struct.
A `StructuralMaterials` is a collection of `Material`s and `Element`s assigning materials to a vector of elements.
### Fields:
- `mats_to_elems` -- Store a dictionary with materials as keys and the corresponding elements as values. """
struct StructuralMaterials{M<:AbstractMaterial,E<:AbstractElement}
    mats_to_elems::Dictionary{M,Vector{E}}
    function StructuralMaterials(mats_to_elems::Dictionary{M,Vector{E}}) where {M<:AbstractMaterial,E<:AbstractElement}
        @assert _element_material_is_unique(mats_to_elems) throw(ArgumentError("Each element must have a single material"))
        new{M,E}(mats_to_elems)
    end
end

"Constructor for empty `StructuralMaterials` with a `Vector` of materials `vmats`. This will assign an empty `Vector` of `Element`s to each material."
StructuralMaterials(vmats::Vector{M}) where {M<:AbstractMaterial} =
    StructuralMaterials(dictionary(map(mat -> mat => Vector{AbstractElement}(), vmats)))

"Returns the `Material` mapped with the label `l`."
Base.getindex(sm::StructuralMaterials, l::L) where {L<:Union{Symbol,String}} =
    first(collect(filter(m -> label(m) == Symbol(l), keys(sm.mats_to_elems))))

"Returns the `Vector` of `Element`s that are conformed by the `Material `m`."
Base.getindex(sm::StructuralMaterials, m::M) where {M<:AbstractMaterial} = sm.mats_to_elems[m]

"Returns the `Vector` of `Material` of the element `e`."
Base.getindex(sm::StructuralMaterials, e::E) where {E<:AbstractElement} =
    first([m for (m, es) in pairs(sm.mats_to_elems) if e in es])

"Returns `Pair`s of `Material` and `Element` in the `StructuralMaterials` `sm`."
Base.pairs(sm::StructuralMaterials) = pairs(sm.mats_to_elems)

"Checks that each `Element` has a single `Material` in the dictionary `mat_dict`."
_element_material_is_unique(mat_dict) = length(unique(values(mat_dict))) == length(values(mat_dict))
