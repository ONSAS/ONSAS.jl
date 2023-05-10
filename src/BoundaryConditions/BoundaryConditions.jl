"""
Module defining the boundary conditions implemented.
Two types of boundary conditions are defined Load (Neumann), and Displacements (Dirichlet). 
Overall, each boundary condition consists of a data type with a label, dofs and values into its fields.
"""
module BoundaryConditions

using Reexport

@reexport import ..Elements: dofs
@reexport import ..Utils: label

export AbstractBoundaryCondition, AbstractNeumannBoundaryCondition,
       AbstractDirichletBoundaryCondition, _apply
export DisplacementBoundaryCondition, FixedDofBoundaryCondition, components
export GlobalLoadBoundaryCondition, LocalPressureBoundaryCondition

""" Abstract supertype for all elements.

An `AbstractBoundaryCondition` object facilitates the process of defining:

    - Dirichlet or displacements boundary conditions.
    - Neumann or load boundary conditions.

**Common methods:**

* [`_apply`](@ref)
* [`Base.values`](@ref)
* [`dofs`](@ref)
* [`label`](@ref)
"""
abstract type AbstractBoundaryCondition end

"Apply the boundary condition `bc` to an`entity`."
function _apply(bc::AbstractBoundaryCondition, entity) end

"Return the degrees of freedom symbol where the boundary condition is imposed"
dofs(bc::AbstractBoundaryCondition) = bc.dofs

"Return the boundary condition label"
label(bc::AbstractBoundaryCondition) = bc.name

"Return the values function imposed to the respective degrees of freedom.
This should be a function of time returning a vector with the same size as the `Node` or `Element` `Dof`s."
Base.values(bc::AbstractBoundaryCondition) = bc.values

#================================#
# Dirichlet boundary conditions  #
#================================#

""" Abstract supertype for all Dirichlet boundary conditions."""
abstract type AbstractDirichletBoundaryCondition <: AbstractBoundaryCondition end

#================================#
# Neumann boundary conditions  #
#================================#
""" Abstract supertype for all Neumann boundary conditions."""
abstract type AbstractNeumannBoundaryCondition <: AbstractBoundaryCondition end

"Abstract functor for a `AbstractLoadBoundaryCondition` that evaluates the load at time `t`."
(lbc::AbstractNeumannBoundaryCondition)(t::Real) = values(lbc)(t)

end # module
