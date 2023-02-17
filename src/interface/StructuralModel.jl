"""
Module defining structure entities interface.
"""
module StructuralModel

using Reexport: @reexport

@reexport using Dictionaries
@reexport using ..Materials
@reexport using ..Elements
@reexport using ..BoundaryConditions
using ..Meshes: AbstractMesh
@reexport using ..Meshes
using ..Utils: row_vector, label

@reexport import ..Elements: dofs, nodes
@reexport import ..Meshes: num_dofs, elements, num_elements, num_nodes

export StructuralMaterials
export StructuralBoundaryConditions, node_bcs, element_bcs, displacement_bcs, load_bcs, fixed_dof_bcs
export Structure, mesh, materials, boundary_conditions, free_dofs

""" Structural materials.
A `StructuralMaterials` is a collection of `Materials` and `Elements` assigning materials to a vector of elements.
### Fields:
- `mats_to_elems` -- Store a dictionary with materials as keys and the corresponding elements as values. """
struct StructuralMaterials{M<:AbstractMaterial,E<:AbstractElement}
    mats_to_elems::Dictionary{M,Vector{E}}
    function StructuralMaterials(mats_to_elems::Dictionary{M,Vector{E}}) where {M<:AbstractMaterial,E<:AbstractElement}
        @assert _element_material_is_unique(mats_to_elems) throw(ArgumentError("Each element must have a single material"))
        new{M,E}(mats_to_elems)
    end
end

Base.getindex(sm::StructuralMaterials, l::L) where {L<:Union{Symbol,String}} =
    collect(filter(m -> label(m) == Symbol(l), keys(sm.mats_to_elems)))[1]

Base.getindex(sm::StructuralMaterials, m::M) where {M<:AbstractMaterial} = sm.mats_to_elems[m]

Base.getindex(sm::StructuralMaterials, e::E) where {E<:AbstractElement} =
    [m for (m, es) in pairs(sm.mats_to_elems) if e in es][1]

Base.pairs(sm::StructuralMaterials) = pairs(sm.mats_to_elems)

"Checks that each element has a single material"
_element_material_is_unique(mat_dict) = length(unique(values(mat_dict))) == length(values(mat_dict))

""" Structural boundary conditions.
A `StructuralBoundaryConditions` is a collection of `BoundaryConditions` defining the boundary conditions of the structure.
### Fields:
- `node_bc` -- Maps each boundary conditions for a vector of nodes. 
- `element_bc` -- Maps each boundary conditions for a vector of elements. 
"""
struct StructuralBoundaryConditions{
    NB<:AbstractBoundaryCondition,LB<:AbstractBoundaryCondition,N<:AbstractNode,E<:AbstractElement
}
    node_bcs::Dictionary{NB,Vector{N}}
    element_bcs::Dictionary{LB,Vector{E}}
end

StructuralBoundaryConditions(node_bcs::Dictionary{BC,Vector{N}}) where {BC<:AbstractBoundaryCondition,N<:AbstractNode} =
    StructuralBoundaryConditions(node_bcs, Dictionary{AbstractBoundaryCondition,Vector{AbstractElement}}())

StructuralBoundaryConditions(element_bcs::Dictionary{BC,Vector{E}}) where {BC<:AbstractBoundaryCondition,E<:AbstractElement} =
    StructuralBoundaryConditions(Dictionary{AbstractBoundaryCondition,Vector{AbstractNode}}(), element_bcs)

"Returns a boundary condition with the label `l`"
function Base.getindex(sb::StructuralBoundaryConditions, l::L) where {L<:Union{Symbol,String}}
    filter(bc -> label(bc) == Symbol(l), vcat(collect(keys(node_bcs(sb))), collect(keys(element_bcs(sb)))))[1]
end

"Returns a collection of nodes and elements that are imposed with the boundary condition `bc`"
function Base.getindex(sb::StructuralBoundaryConditions{NB,LB}, bc::BC) where
{NB<:AbstractBoundaryCondition,LB<:AbstractBoundaryCondition,BC<:AbstractBoundaryCondition}

    bc_elements = Vector{Union{AbstractElement,AbstractNode}}()

    BC <: NB && bc ∈ keys(node_bcs(sb)) && push!(bc_elements, node_bcs(sb)[bc]...)
    BC <: LB && bc ∈ keys(element_bcs(sb)) && push!(bc_elements, element_bcs(sb)[bc]...)

    isempty(bc_elements) ? throw(KeyError("Boundary condition $bc not found")) : return bc_elements
end

"Returns a collection of boundary conditions applied to the node"
Base.getindex(sb::StructuralBoundaryConditions, n::AbstractNode) = keys(filter(x -> n ∈ x, node_bcs(sb)))

"Returns a collection of boundary conditions applied to the element"
Base.getindex(sb::StructuralBoundaryConditions, e::AbstractElement) = keys(filter(x -> e ∈ x, element_bcs(sb)))

"Returns boundary conditions applied to nodes"
node_bcs(se::StructuralBoundaryConditions) = se.node_bcs

"Returns load boundary conditions applied to elements"
element_bcs(se::StructuralBoundaryConditions) = se.element_bcs

"Returns displacements boundary conditions"
function displacement_bcs(se::StructuralBoundaryConditions)
    vbc = Vector{DisplacementBoundaryCondition}()
    disp_bc_nodes = filter(bc -> bc isa DisplacementBoundaryCondition, keys(node_bcs(se)))
    push!(vbc, disp_bc_nodes...)
    disp_bc_elements = filter(bc -> bc isa DisplacementBoundaryCondition, keys(element_bcs(se)))
    push!(vbc, disp_bc_elements...)
    return unique(vbc)
end

"Returns fixed dof boundary conditions"
function fixed_dof_bcs(se::StructuralBoundaryConditions)
    vbc = Vector{FixedDofBoundaryCondition}()
    fixed_bc_nodes = filter(bc -> bc isa FixedDofBoundaryCondition, keys(node_bcs(se)))
    push!(vbc, fixed_bc_nodes...)
    fixed_bc_elements = filter(bc -> bc isa FixedDofBoundaryCondition, keys(element_bcs(se)))
    push!(vbc, fixed_bc_elements...)
    return unique(vbc)
end

"Returns load boundary conditions"
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
* [`free_dof_indexes`](@ref)
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
"Returns the mesh of the structure `s`"
mesh(s::AbstractStructure) = s.mesh

"Returns structural mesh dofs of the structure `s`"
dofs(s::AbstractStructure) = dofs(mesh(s))

"Returns number of dofs of the structure `s`"
num_dofs(s::AbstractStructure) = num_dofs(mesh(s))

"Returns free dofs (or not fixed) dofs of the structure"
free_dofs(s::AbstractStructure) = s.free_dofs

"Returns structural mesh nodes of the structure `s`"
nodes(s::AbstractStructure) = nodes(mesh(s))

"Returns the number of nodes of the structure `s`"
num_nodes(s::AbstractStructure) = num_nodes(mesh(s))

"Returns structural mesh elements"
elements(s::AbstractStructure) = elements(mesh(s))

"Returns number of elements of the structure `s`"
num_elements(s::AbstractStructure) = num_elements(mesh(s))

# Boundary Conditions 
"Returns the structural boundary conditions of the structure `s`"
boundary_conditions(s::AbstractStructure) = s.bcs

"Returns displacement boundary conditions"
displacement_bcs(s::AbstractStructure) = displacement_bcs(s.bcs)

"Returns load boundary conditions"
load_bcs(s::AbstractStructure) = load_bcs(s.bcs)

"Returns nodal boundary conditions"
node_bcs(s::AbstractStructure) = node_bcs(s.bcs)

"Returns element boundary conditions"
element_bcs(s::AbstractStructure) = element_bcs(s.bcs)

# Materials
"Returns the materials of the structure `s`"
materials(s::AbstractStructure) = s.materials

# Structure 
"""
An `Structure` object facilitates the process of assembling and creating the structural analysis. 
### Fields:
- `mesh`             -- Stores the structure mesh. 
- `materials`        -- Stores the structural materials of the structure. 
- `elements`         -- Stores the structural elements of the structure.
- `bcs`              -- Stores the structural boundary conditions of the structure.
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

function Structure(
    mesh::AbstractMesh{dim},
    materials::StructuralMaterials{M,E},
    bcs::StructuralBoundaryConditions{NB,LB},
) where {dim,M,E,NB,LB}

    default_free_dofs = Vector{Dof}()
    for node_dofs in dofs(mesh)
        [push!(default_free_dofs, vec_dof...) for vec_dof in collect(values(node_dofs))]
    end
    return Structure(mesh, materials, bcs, default_free_dofs)
end

end # module
