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
using ..Utils

@reexport import ..Meshes: dofs, num_dofs, nodes, num_nodes, faces, num_faces,
                           elements, num_elements
@reexport import ..Handlers: mesh
export Structure, materials, boundary_conditions, num_free_dofs, free_dofs

# Structural properties 
include("./StructuralMaterials.jl")
include("./StructuralBoundaryConditions.jl")
include("./StructuralEntities.jl")

# ==========
# Structure
# ==========

""" 
Abstract supertype to define a new structure. 
An `AbstractStructure` object facilitates the process of assigning materials, elements, 
and boundary conditions to the mesh. Moreover this struct is used to solve an specific 
structural analysis. 

**Common methods:**

## Mesh: 
* [`dofs`](@ref)
* [`num_dofs`](@ref)
* [`free_dofs`](@ref)
* [`num_free_dofs`](@ref)
* [`nodes`](@ref)
* [`num_nodes`](@ref)
* [`faces`](@ref)
* [`num_faces`](@ref)
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
"Return the `Mesh` of the `AbstractStructure` `s`"
mesh(s::AbstractStructure) = s.mesh

"Return the `Dof`s of the `AbstractStructure` `s`"
dofs(s::AbstractStructure) = dofs(mesh(s))

"Return the number of `Dof`s of the `AbstractStructure` `s`"
num_dofs(s::AbstractStructure) = num_dofs(mesh(s))

"Return free `Dof`s of the structure `s`"
free_dofs(s::AbstractStructure) = s.free_dofs

"Return the number of free `Dof`s of the `AbstractStructure` `s`"
num_free_dofs(s::AbstractStructure) = length(free_dofs(s))

"Return a `Vector` of `Node`s defined in the `AbstractStructure` `s`."
nodes(s::AbstractStructure) = nodes(mesh(s))

"Return the number of `Node`s of the `AbstractStructure` `s`"
num_nodes(s::AbstractStructure) = num_nodes(mesh(s))

"Return the `Element`s of the `AbstractStructure` `s`"
elements(s::AbstractStructure) = elements(mesh(s))

"Return the number of `Element`s of the `AbstractStructure` `s`"
num_elements(s::AbstractStructure) = num_elements(mesh(s))

"Return the `Face`s of the `AbstractStructure` `s`"
faces(s::AbstractStructure) = faces(mesh(s))

"Return the number of `Face`s of the `AbstractStructure` `s`"
num_faces(s::AbstractStructure) = num_faces(mesh(s))

# Boundary Conditions 
"Return the `StructuralBoundaryConditions` of the `AbstractStructure` `s`"
boundary_conditions(s::AbstractStructure) = s.bcs

"Return the `DisplacementBoundaryCondition`s of the `AbstractStructure` `s`"
displacement_bcs(s::AbstractStructure) = displacement_bcs(s.bcs)

"Return the `LocalLoad`s of the `AbstractStructure` `s`"
load_bcs(s::AbstractStructure) = load_bcs(s.bcs)

"Return the `AbstractBoundaryCondition`s imposed to `Node`s in the `AbstractStructure` `s`"
node_bcs(s::AbstractStructure) = node_bcs(s.bcs)

"Return the `AbstractBoundaryCondition`s imposed to `Element`s in the `AbstractStructure` `s`"
element_bcs(s::AbstractStructure) = element_bcs(s.bcs)

# Materials
"Return the `StructuralMaterials`s of the `AbstractStructure` `s`"
materials(s::AbstractStructure) = s.materials

#====================================#
# Abstract structure implementations #
#====================================#
include("./Structure.jl")

end # module
