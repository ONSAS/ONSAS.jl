""" Structural boundary conditions.
A `StructuralBoundaryConditions` is a collection of `BoundaryConditions` defining the boundary conditions of the structure.
### Fields:
- `node_bc`    -- Maps each boundary conditions for a vector of nodes. 
- `element_bc` -- Maps each boundary conditions for a vector of elements. 
"""
struct StructuralBoundaryConditions{
    NB<:AbstractBoundaryCondition,NE<:AbstractBoundaryCondition,N<:AbstractNode,E<:AbstractElement
}
    node_bcs::Dictionary{NB,Vector{N}}
    element_bcs::Dictionary{NE,Vector{E}}
end

"Constructor for `StructuralBoundaryConditions` with node boundary conditions."
StructuralBoundaryConditions(node_bcs::Dictionary{BC,Vector{N}}) where {BC<:AbstractBoundaryCondition,N<:AbstractNode} =
    StructuralBoundaryConditions(node_bcs, Dictionary{AbstractBoundaryCondition,Vector{AbstractElement}}())

"Constructor for `StructuralBoundaryConditions` with element boundary conditions."
StructuralBoundaryConditions(element_bcs::Dictionary{BC,Vector{E}}) where {BC<:AbstractBoundaryCondition,E<:AbstractElement} =
    StructuralBoundaryConditions(Dictionary{AbstractBoundaryCondition,Vector{AbstractNode}}(), element_bcs)

"Returns the `BoundaryCondition` with the label `l` in the `StructuralBoundaryConditions` `sb`."
function Base.getindex(sb::StructuralBoundaryConditions, l::L) where {L<:Union{Symbol,String}}
    filter(bc -> label(bc) == Symbol(l), vcat(collect(keys(node_bcs(sb))), collect(keys(element_bcs(sb)))))[1]
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

"Returns the `Dictionary` of `BoundaryConditions`s applied to `Element`s."
element_bcs(se::StructuralBoundaryConditions) = se.element_bcs

"Returns a `Vector` of `DisplacementBoundaryCondition`s applied to `Node`s and `Element`s in the `StructuralBoundaryConditions` `se`."
function displacement_bcs(se::StructuralBoundaryConditions)
    vbc = Vector{DisplacementBoundaryCondition}()
    disp_bc_nodes = filter(bc -> bc isa DisplacementBoundaryCondition, keys(node_bcs(se)))
    push!(vbc, disp_bc_nodes...)
    disp_bc_elements = filter(bc -> bc isa DisplacementBoundaryCondition, keys(element_bcs(se)))
    push!(vbc, disp_bc_elements...)
    return unique(vbc)
end

"Returns a `Vector` of `FixedDofBoundaryCondition`s applied to `Node`s and `Element`s. in the `StructuralBoundaryConditions` `se`."
function fixed_dof_bcs(se::StructuralBoundaryConditions)
    vbc = Vector{FixedDofBoundaryCondition}()
    fixed_bc_nodes = filter(bc -> bc isa FixedDofBoundaryCondition, keys(node_bcs(se)))
    push!(vbc, fixed_bc_nodes...)
    fixed_bc_elements = filter(bc -> bc isa FixedDofBoundaryCondition, keys(element_bcs(se)))
    push!(vbc, fixed_bc_elements...)
    return unique(vbc)
end

"Returns a `Vector` of `FixedDofBoundaryCondition`s applied to `Node`s and `Element`s in the `StructuralBoundaryConditions` `se`."
function load_bcs(se::StructuralBoundaryConditions)
    vbc = Vector{AbstractLoadBoundaryCondition}()
    load_bc_nodes = filter(bc -> bc isa AbstractLoadBoundaryCondition, keys(node_bcs(se)))
    push!(vbc, load_bc_nodes...)
    load_bc_elements = filter(bc -> bc isa AbstractLoadBoundaryCondition, keys(element_bcs(se)))
    push!(vbc, load_bc_elements...)
    return unique(vbc)
end