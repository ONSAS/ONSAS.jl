using Reexport
using Dictionaries

using ..Materials
using ..Elements
using ..Utils

@reexport import ..Elements: apply!

export StructuralMaterials, element_materials

""" 
A `StructuralMaterials` is a collection of `Material`s and `Element`s assigning materials to a vector of elements.
"""
struct StructuralMaterials{M<:AbstractMaterial,E<:AbstractElement}
    "Stores materials as keys and the corresponding elements as values"
    mats_to_elems::Dictionary{M,Vector{E}}
    function StructuralMaterials(mats_to_elems::Dictionary{M,Vector{E}}) where
             {M<:AbstractMaterial,E<:AbstractElement}
        @assert _element_material_is_unique(mats_to_elems) error("Each element must have a single material")
        # Abstract is used to replace materials 
        new{AbstractMaterial,E}(mats_to_elems)
    end
end

"Return a `Dictionary` with `Material`s as keys and the corresponding `Element`s as values."
element_materials(sm::StructuralMaterials) = sm.mats_to_elems

"Constructor for empty `StructuralMaterials` with a `Vector` of materials `vmats`. This will assign an empty `Vector` of `Element`s to each material."
function StructuralMaterials(vmats::Vector{M}) where {M<:AbstractMaterial}
    StructuralMaterials(dictionary(map(mat -> mat => Vector{AbstractElement}(), vmats)))
end

StructuralMaterials(vmats::AbstractMaterial...) = StructuralMaterials(collect(vmats))

"Return the `Material` mapped with the label `l`."
function Base.getindex(sm::StructuralMaterials, l::L) where {L<:Union{Symbol,AbstractString}}
    materials_label_l = collect(filter(m -> label(m) == Symbol(l), keys(element_materials(sm))))
    @assert length(materials_label_l) == 1 throw(ArgumentError("The label $l is not found."))
    first(materials_label_l)
end

"Return the `Vector` of `Element`s that are conformed by the `Material `m`."
Base.getindex(sm::StructuralMaterials, m::M) where {M<:AbstractMaterial} = element_materials(sm)[m]

"Return the `Vector` of `Material` of the element `e`."
function Base.getindex(sm::StructuralMaterials, e::E) where {E<:AbstractElement}
    first([m for (m, es) in pairs(element_materials(sm)) if e in es])
end

"Return `Pair`s of `Material` and `Element` in the `StructuralMaterials` `sm`."
Base.pairs(sm::StructuralMaterials) = pairs(element_materials(sm))

"Checks that each `Element` has a single `Material` in the dictionary `mat_dict`."
_element_material_is_unique(mat_dict) = length(unique(values(mat_dict))) == length(values(mat_dict))

"Apply the `StructuralBoundaryConditions` to the `AbstractMesh` `m`. For this is required sets 
into the `Mesh` and the corresponding boundary condition labels declared in `bcs`."
function apply!(sm::StructuralMaterials, m::AbstractMesh)
    vec_elements = elements(m)
    element_sets = element_set(m)
    for (mat, elements) in pairs(element_materials(sm))
        mat_label = string(label(mat))
        [push!(elements, vec_elements[element_index]) for element_index in element_sets[mat_label]]
    end
end

"Delete the `Material` `m` from the `StructuralMaterials` `sm`."
function Base.delete!(sm::StructuralMaterials, m::M) where {M<:AbstractMaterial}
    delete!(element_materials(sm), m)
end

"Insert the `Material` `m` into the `StructuralMaterials` `sm`."
function Base.insert!(sm::StructuralMaterials, m::M,
                      mat_elements::Vector{<:AbstractElement}=Vector{AbstractElement}()) where {M<:AbstractMaterial}
    insert!(element_materials(sm), m, mat_elements)
end

"Replace the `AbstractMaterial` with the label `l` for the new material `new_m` in the `StructuralMaterials` `sm`.
The previous material `Element`s are assigned to the new."
function Base.replace!(sm::StructuralMaterials,
                       new_material::AbstractMaterial,
                       label::Label=label(new_material))
    old_material = sm[label]
    material_elements = sm[old_material]
    delete!(sm, old_material)
    insert!(sm, new_material, material_elements)
end
