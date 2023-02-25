using ..Elements: AbstractNode, AbstractFace, AbstractElement
using Reexport: @reexport

@reexport import ..BoundaryConditions: apply

export StructuralBoundaryConditions, all_bcs, node_bcs, face_bcs, element_bcs, displacement_bcs, load_bcs, fixed_dof_bcs

""" Structural boundary conditions.
A `StructuralBoundaryConditions` is a collection of `BoundaryConditions` defining the boundary conditions of the structure.
### Fields:
- `node_bcs`    -- Maps each boundary conditions for a vector of nodes. 
- `face_bcs`    -- Maps each boundary conditions for a vector of faces. 
- `element_bcs` -- Maps each boundary conditions for a vector of elements. 
"""
Base.@kwdef struct StructuralBoundaryConditions{
    NB<:AbstractBoundaryCondition,FB<:AbstractBoundaryCondition,EB<:AbstractBoundaryCondition,
    N<:AbstractNode,F<:AbstractFace,E<:AbstractElement
}
    node_bcs::Dictionary{NB,Vector{N}} = Dictionary{AbstractBoundaryCondition,Vector{AbstractNode}}()
    face_bcs::Dictionary{FB,Vector{F}} = Dictionary{AbstractBoundaryCondition,Vector{AbstractFace}}()
    element_bcs::Dictionary{EB,Vector{E}} = Dictionary{AbstractBoundaryCondition,Vector{AbstractElement}}()
end

"Returns the `BoundaryCondition` with the label `l` in the `StructuralBoundaryConditions` `sb`."
function Base.getindex(sb::StructuralBoundaryConditions, l::L) where {L<:Union{Symbol,String}}
    first(filter(bc -> label(bc) == Symbol(l), all_bcs(sb)))
end

"Returns the `Vector` of `Node`s and `Element`s where the `BoundaryCondition` `bc` is imposed."
function Base.getindex(sb::StructuralBoundaryConditions{NB,NF,EB}, bc::BC) where
{NB<:AbstractBoundaryCondition,NF<:AbstractBoundaryCondition,EB<:AbstractBoundaryCondition,BC<:AbstractBoundaryCondition}

    bc_entities = Vector{Union{AbstractElement,AbstractNode,AbstractFace}}()

    BC <: NB && bc ∈ keys(node_bcs(sb)) && push!(bc_entities, node_bcs(sb)[bc]...)
    BC <: NF && bc ∈ keys(face_bcs(sb)) && push!(bc_entities, face_bcs(sb)[bc]...)
    BC <: EB && bc ∈ keys(element_bcs(sb)) && push!(bc_entities, element_bcs(sb)[bc]...)

    isempty(bc_entities) ? throw(KeyError("Boundary condition $bc not found")) : return bc_entities
end

"Returns the `Vector` of `BoundaryConditions`s applied to the `AbstractNode` `n`."
Base.getindex(sb::StructuralBoundaryConditions, n::AbstractNode) = keys(filter(x -> n ∈ x, node_bcs(sb)))

"Returns the `Vector` of `BoundaryConditions`s applied to the `AbstractFace` `f`."
Base.getindex(sb::StructuralBoundaryConditions, f::AbstractFace) = keys(filter(x -> f ∈ x, face_bcs(sb)))

"Returns the `Vector` of `BoundaryConditions`s applied to the `AbstractElement` `e`."
Base.getindex(sb::StructuralBoundaryConditions, e::AbstractElement) = keys(filter(x -> e ∈ x, element_bcs(sb)))

"Returns the dictionary of `BoundaryConditions`s applied to `Node`s."
node_bcs(se::StructuralBoundaryConditions) = se.node_bcs

"Returns the dictionary of `BoundaryConditions`s applied to `Face`s."
face_bcs(se::StructuralBoundaryConditions) = se.face_bcs

"Returns the `Dictionary` of `BoundaryConditions`s applied to `Element`s."
element_bcs(se::StructuralBoundaryConditions) = se.element_bcs

"Returns all `BoundaryConditions`s defined into `StructuralBoundaryConditions`."
all_bcs(se::StructuralBoundaryConditions) = unique(
    vcat(collect(keys(node_bcs(se))), collect(keys(face_bcs(se))), collect(keys(element_bcs(se))))
)

"Returns a `Vector` of `DisplacementBoundaryCondition`s applied to `Node`s and `Element`s in the `StructuralBoundaryConditions` `se`."
function displacement_bcs(se::StructuralBoundaryConditions)
    vbc = Vector{DisplacementBoundaryCondition}()
    disp_bc = filter(bc -> bc isa DisplacementBoundaryCondition, all_bcs(se))
    push!(vbc, disp_bc...)
    return unique(vbc)
end

"Returns a `Vector` of `FixedDofBoundaryCondition`s applied to `Node`s, `Face`s and `Element`s. in the `StructuralBoundaryConditions` `se`."
function fixed_dof_bcs(se::StructuralBoundaryConditions)
    vbc = Vector{FixedDofBoundaryCondition}()
    fixed_bcs = filter(bc -> bc isa FixedDofBoundaryCondition, all_bcs(se))
    push!(vbc, fixed_bcs...)
    return unique(vbc)
end

"Returns a `Vector` of `Dof`s to delete given a `FixedDofBoundaryCondition` and a set of `StructuralBoundaryConditions` `bcs`."
function apply(bcs::StructuralBoundaryConditions, fbc::FixedDofBoundaryCondition)
    # Extract nodes, faces and elements 
    entities = bcs[fbc]
    dofs_to_delete = Dof[]
    [push!(dofs_to_delete, apply(fbc, e)...) for e in entities]
    return unique!(dofs_to_delete)
end

"Returns a `Vector` of `FixedDofBoundaryCondition` `f_bcs` and a set of `StructuralBoundaryConditions` `bcs`."
function apply(bcs::StructuralBoundaryConditions, f_bcs::Vector{<:FixedDofBoundaryCondition})
    dofs_to_delete = Dof[]
    [push!(dofs_to_delete, compute_fixed_dofs(bcs, fbc)...) for fbc in f_bcs]
    return unique!(dofs_to_delete)
end

"Returns a `Vector` of `FixedDofBoundaryCondition`s applied to `Node`s and `Element`s in the `StructuralBoundaryConditions` `bcs`."
function load_bcs(bcs::StructuralBoundaryConditions)
    vbc = Vector{AbstractLoadBoundaryCondition}()
    load_bcs = filter(bc -> bc isa AbstractLoadBoundaryCondition, all_bcs(bcs))
    push!(vbc, load_bcs...)
    return unique(vbc)
end

"Returns a `Vector` of `Dof`s and values to apply a `LoadDofBoundaryCondition` and a set of `StructuralBoundaryConditions` `bcs`
at time `t`."
function apply(bcs::StructuralBoundaryConditions, lbc::AbstractLoadBoundaryCondition, t::Real)
    # Extract nodes, faces and elements 
    entities = bcs[lbc]
    dofs_to_load = Dof[]
    load_vec = Float64[]
    for e in entities
        dofs_lbc_e, load_vec_e = apply(lbc, e, t)
        push!(load_vec, load_vec_e...)
        push!(dofs_to_load, dofs_lbc_e...)
    end
    return dofs_to_load, load_vec
end