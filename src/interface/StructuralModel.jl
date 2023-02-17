"""
Module defining structure entities interface.
"""
module StructuralModel

using Reexport: @reexport

@reexport using Dictionaries: Dictionary
@reexport using ..Materials
@reexport using ..Elements
@reexport using ..BoundaryConditions
@reexport using ..Meshes
using ..Utils: row_vector, label

import ..Elements: dofs, nodes
import ..Meshes: num_dofs, elements, num_elements, num_nodes

export StructuralMaterials
export StructuralBoundaryConditions, node_bcs, element_bcs, displacement_bcs, load_bcs, fixed_dof_bcs
export Structure, mesh, materials, boundary_conditions, num_free_dofs, free_dofs

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

"Returns the `Material` mapped with the label `l`."
Base.getindex(sm::StructuralMaterials, l::L) where {L<:Union{Symbol,String}} =
    collect(filter(m -> label(m) == Symbol(l), keys(sm.mats_to_elems)))[1]

"Returns the `Vector` of `Element`s that are conformed by the `Material `m`."
Base.getindex(sm::StructuralMaterials, m::M) where {M<:AbstractMaterial} = sm.mats_to_elems[m]

"Returns the `Vector` of `Material` of the element `e`."
Base.getindex(sm::StructuralMaterials, e::E) where {E<:AbstractElement} =
    [m for (m, es) in pairs(sm.mats_to_elems) if e in es][1]

"Returns `Pair`s of `Material` and `Element` in the `StructuralMaterials` `sm`."
Base.pairs(sm::StructuralMaterials) = pairs(sm.mats_to_elems)

"Checks that each `Element` has a single `Material` in the dictionary `mat_dict`."
_element_material_is_unique(mat_dict) = length(unique(values(mat_dict))) == length(values(mat_dict))

""" Structural boundary conditions.
A `StructuralBoundaryConditions` is a collection of `BoundaryConditions` defining the boundary conditions of the structure.
### Fields:
- `node_bc`    -- Maps each boundary conditions for a vector of nodes. 
- `element_bc` -- Maps each boundary conditions for a vector of elements. 
"""
struct StructuralBoundaryConditions{
    NB<:AbstractBoundaryCondition,LB<:AbstractBoundaryCondition,N<:AbstractNode,E<:AbstractElement
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

# ==========
# Structure
# ==========

""" Abstract supertype to define a new structure. An `AbstractStructure` object facilitates the process of assigning materials, 
elements, and boundary conditions to the mesh. Moreover this struct is used to solve an specific structural analysis. 
**Common methods:**

## Mesh: 
* [`dofs`](@ref)
* [`num_dofs`](@ref)
* [`free_dofs`](@ref)
* [`num_free_dofs`](@ref)
* [`nodes`](@ref)
* [`num_nodes`](@ref)
* [`elements`](@ref)
* [`num_elements`](@ref)
* [`mesh`](@ref)

## Boundary conditions:
* [`boundary_conditions`](@ref)
* [`displacement_bcs`](@ref)
* [`load_bcs`](@ref)
* [`node_bcs`](@ref)
* [`element_bcs`](@ref)

## Materials:   
* [`materials`](@ref)

"""
abstract type AbstractStructure{dim,M,E} end

# Mesh methods
"Returns the `Mesh` of the `AbstractStructure` `s`"
mesh(s::AbstractStructure) = s.mesh

"Returns the `Dof`s of the `AbstractStructure` `s`"
dofs(s::AbstractStructure) = dofs(mesh(s))

"Returns the number of `Dof`s of the `AbstractStructure` `s`"
num_dofs(s::AbstractStructure) = num_dofs(mesh(s))

"Returns free `Dof`s of the structure `s`"
free_dofs(s::AbstractStructure) = s.free_dofs

"Returns the number of free `Dof`s of the `AbstractStructure` `s`"
num_free_dofs(s::AbstractStructure) = length(free_dofs(s))

"Returns the `Node`s of the `AbstractStructure` `s`"
nodes(s::AbstractStructure) = nodes(mesh(s))

"Returns the number of `Node`s of the `AbstractStructure` `s`"
num_nodes(s::AbstractStructure) = num_nodes(mesh(s))

"Returns the `Element`s of the `AbstractStructure` `s`"
elements(s::AbstractStructure) = elements(mesh(s))

"Returns the number of `Element`s of the `AbstractStructure` `s`"
num_elements(s::AbstractStructure) = num_elements(mesh(s))

# Boundary Conditions 
"Returns the `StructuralBoundaryConditions` of the `AbstractStructure` `s`"
boundary_conditions(s::AbstractStructure) = s.bcs

"Returns the `DisplacementBoundaryCondition`s of the `AbstractStructure` `s`"
displacement_bcs(s::AbstractStructure) = displacement_bcs(s.bcs)

"Returns the `LocalLoadBoundaryCondition`s of the `AbstractStructure` `s`"
load_bcs(s::AbstractStructure) = load_bcs(s.bcs)

"Returns the `AbstractBoundaryCondition`s imposed to `Node`s in the `AbstractStructure` `s`"
node_bcs(s::AbstractStructure) = node_bcs(s.bcs)

"Returns the `AbstractBoundaryCondition`s imposed to `Element`s in the `AbstractStructure` `s`"
element_bcs(s::AbstractStructure) = element_bcs(s.bcs)

# Materials
"Returns the `StructuralMaterials`s of the `AbstractStructure` `s`"
materials(s::AbstractStructure) = s.materials

# Structure 
"""
An `Structure` object facilitates the process of assembling and creating the structural analysis. 
### Fields:
- `mesh`      -- Stores the structural mesh. 
- `materials` -- Stores the structural materials of the structure. 
- `elements`  -- Stores the structural elements of the structure.
- `bcs`       -- Stores the structural boundary conditions of the structure.
- `free_dofs` -- Stores the free degrees of freedom.
"""
mutable struct Structure{dim,M,E,NB,LB} <: AbstractStructure{dim,M,E}
    mesh::AbstractMesh{dim}
    materials::StructuralMaterials{M,E}
    bcs::StructuralBoundaryConditions{NB,LB}
    free_dofs::Vector{Dof}
    function Structure(
        mesh::AbstractMesh{dim},
        materials::StructuralMaterials{M,E},
        bcs::StructuralBoundaryConditions{NB,LB},
        free_dofs::Vector{Dof}
    ) where {dim,M,E,B,NB,LB}
        return new{dim,M,E,NB,LB}(mesh, materials, bcs, free_dofs)
    end
end

"Constructor with  `StructuralMaterials` `materials`,  `StructuralBoundaryConditions` `bcs` 
and `AbstractMesh` `mesh` seting fixed dofs with `FixedDofBoundaryCondition` defined in `bcs`"
function Structure(
    mesh::AbstractMesh{dim},
    materials::StructuralMaterials{M,E},
    bcs::StructuralBoundaryConditions{NB,LB},
) where {dim,M,E,NB,LB}

    default_free_dofs = Vector{Dof}()
    for node_dofs in dofs(mesh)
        [push!(default_free_dofs, vec_dof...) for vec_dof in collect(values(node_dofs))]
    end

    fixed_dofs = compte_fixed_dofs(bcs, fixed_dof_bcs(bcs))

    deleteat!(default_free_dofs, findall(x -> x in fixed_dofs, default_free_dofs))

    return Structure(mesh, materials, bcs, default_free_dofs)
end

#### BoundaryConditions Applied to the structure
"Computes `Dof`s to delete given a `FixedDofBoundaryCondition` and a set of `StructuralBoundaryConditions` `bcs`."
function compte_fixed_dofs(bcs::StructuralBoundaryConditions, fbc::FixedDofBoundaryCondition)

    # Extract dofs to apply the bc
    fbc_dofs_symbols = dofs(fbc)

    # Extract nodes and elements 
    entities = bcs[fbc]

    dofs_to_delete = Dof[]

    for dof_symbol in fbc_dofs_symbols
        dofs_entities = getindex.(dofs.(entities), dof_symbol)
        for component in components(fbc)
            dofs_to_delete_fbc = getindex.(dofs_entities, component)
            push!(dofs_to_delete, dofs_to_delete_fbc...)
        end
    end
    return dofs_to_delete
end

"Applies a `Vector` of `FixedDofBoundaryCondition` `f_bcs` and a set of `StructuralBoundaryConditions` `bcs`."
function compte_fixed_dofs(bcs::StructuralBoundaryConditions, f_bcs::Vector{<:FixedDofBoundaryCondition})
    dofs_to_delete = Dof[]
    [push!(dofs_to_delete, compte_fixed_dofs(bcs, fbc)...) for fbc in f_bcs]
    return dofs_to_delete
end

end # module
