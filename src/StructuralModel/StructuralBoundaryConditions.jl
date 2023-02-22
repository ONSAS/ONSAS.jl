using ..Elements: AbstractNode, AbstractFace, AbstractElement

export StructuralBoundaryConditions, all_bcs, node_bcs, face_bcs, element_bcs, displacement_bcs, load_bcs, fixed_dof_bcs

""" Structural boundary conditions.
A `StructuralBoundaryConditions` is a collection of `BoundaryConditions` defining the boundary conditions of the structure.
### Fields:
- `node_bcs`    -- Maps each boundary conditions for a vector of nodes. 
- `face_bcs`    -- Maps each boundary conditions for a vector of faces. 
- `element_bcs` -- Maps each boundary conditions for a vector of elements. 
"""
Base.@kwdef struct StructuralBoundaryConditions{
    NB<:AbstractBoundaryCondition,NF<:AbstractBoundaryCondition,NE<:AbstractBoundaryCondition,
    N<:AbstractNode,F<:AbstractFace,E<:AbstractElement
}
    node_bcs::Dictionary{NB,Vector{N}} = Dictionary{AbstractBoundaryCondition,Vector{AbstractNode}}()
    face_bcs::Dictionary{NF,Vector{F}} = Dictionary{AbstractBoundaryCondition,Vector{AbstractFace}}()
    element_bcs::Dictionary{NE,Vector{E}} = Dictionary{AbstractBoundaryCondition,Vector{AbstractElement}}()
end

"Returns the `BoundaryCondition` with the label `l` in the `StructuralBoundaryConditions` `sb`."
function Base.getindex(sb::StructuralBoundaryConditions, l::L) where {L<:Union{Symbol,String}}
    filter(bc -> label(bc) == Symbol(l), vcat(
        collect(keys(node_bcs(sb))), collect(keys(face_bcs(sb))), collect(keys(element_bcs(sb)))))[1]
end

"Returns the `Vector` of `Node`s and `Element`s where the `BoundaryCondition` `bc` is imposed."
function Base.getindex(sb::StructuralBoundaryConditions{NB,LB}, bc::BC) where
{NB<:AbstractBoundaryCondition,LB<:AbstractBoundaryCondition,BC<:AbstractBoundaryCondition}

    bc_elements = Vector{Union{AbstractElement,AbstractNode}}()

    BC <: NB && bc ∈ keys(node_bcs(sb)) && push!(bc_elements, node_bcs(sb)[bc]...)
    BC <: LB && bc ∈ keys(element_bcs(sb)) && push!(bc_elements, element_bcs(sb)[bc]...)

    isempty(bc_elements) ? throw(KeyError("Boundary condition $bc not found")) : return bc_elements
end

"Returns the `Vector` of `BoundaryConditions`s applied to the `Node` `n`."
Base.getindex(sb::StructuralBoundaryConditions, n::AbstractNode) = keys(filter(x -> n ∈ x, node_bcs(sb)))

"Returns the `Vector` of `BoundaryConditions`s applied to the `Element` `e`."
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

"Returns a `Vector` of `FixedDofBoundaryCondition`s applied to `Node`s and `Element`s in the `StructuralBoundaryConditions` `se`."
function load_bcs(se::StructuralBoundaryConditions)
    vbc = Vector{AbstractLoadBoundaryCondition}()
    load_bcs = filter(bc -> bc isa AbstractLoadBoundaryCondition, all_bcs(se))
    push!(vbc, load_bcs...)
    return unique(vbc)
end