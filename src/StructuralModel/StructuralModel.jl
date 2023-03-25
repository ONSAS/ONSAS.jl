"""
Module defining structure entities interface.
Each Structure consists of a data type with materials, boundary conditions, mesh and a vector of free dofs.
"""
module StructuralModel

using Reexport: @reexport

@reexport using ..Materials
@reexport using ..Elements
@reexport using ..BoundaryConditions
@reexport using ..Meshes
using ..Utils: row_vector, label

import ..Elements: dofs, nodes
import ..Meshes: mesh, num_dofs, elements, num_elements, num_nodes

export Structure, materials, boundary_conditions, num_free_dofs, free_dofs


# Structural properties 
include("./StructuralMaterials.jl")
include("./StructuralBoundaryConditions.jl")
include("./StructuralEntities.jl")


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

"Returns a `Vector` of `Node`s defined in the `AbstractStructure` `s`."
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

"Returns the `LocalPressureBoundaryCondition`s of the `AbstractStructure` `s`"
load_bcs(s::AbstractStructure) = load_bcs(s.bcs)

"Returns the `AbstractBoundaryCondition`s imposed to `Node`s in the `AbstractStructure` `s`"
node_bcs(s::AbstractStructure) = node_bcs(s.bcs)

"Returns the `AbstractBoundaryCondition`s imposed to `Element`s in the `AbstractStructure` `s`"
element_bcs(s::AbstractStructure) = element_bcs(s.bcs)

# Materials
"Returns the `StructuralMaterials`s of the `AbstractStructure` `s`"
materials(s::AbstractStructure) = s.materials

#====================================#
# Abstract structure implementations #
#====================================#
include("./Structure.jl")


end # module
