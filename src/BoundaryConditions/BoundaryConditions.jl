"""
Module defining the boundary conditions implemented.
Two types of boundary conditions are defined Load (Neumann), and Displacements (Dirichlet). 
Overall, each boundary condition consists of a data type with a label, dofs and values into its fields.
"""
module BoundaryConditions

using Reexport: @reexport

@reexport import ..Elements: dofs
@reexport import ..Utils: label

export AbstractBoundaryCondition, AbstractDisplacementBoundaryCondition, AbstractLoadBoundaryCondition, _apply
export DisplacementBoundaryCondition, FixedDofBoundaryCondition, components
export GlobalLoadBoundaryCondition, LocalLoadBoundaryCondition


""" Abstract supertype for all elements.

An `AbstractBoundaryCondition` object facilitates the process of defining:

    - Displacements (Dirichlet) boundary conditions.
    - Load (Neumann) boundary conditions.

**Common methods:**

* [`dofs`](@ref)
* [`label`](@ref)
"""
abstract type AbstractBoundaryCondition end

"Applies the boundary condition `bc` to an`entity`."
function _apply(bc::AbstractBoundaryCondition, entity) end

"Returns the degrees of freedom symbol where the boundary condition is imposed"
dofs(bc::AbstractBoundaryCondition) = bc.dofs

"Returns the boundary condition label"
label(bc::AbstractBoundaryCondition) = bc.name

"Returns the values function imposed to the respective degrees of freedom.
This should be a function of time returning a vector with the same size as the node or element dofs."
Base.values(bc::AbstractBoundaryCondition) = bc.values

# ================================
# Displacement Boundary Conditions 
# ================================

""" Abstract supertype for all displacement boundary conditions."""
abstract type AbstractDisplacementBoundaryCondition <: AbstractBoundaryCondition end

include("FixedDofBoundaryCondition.jl")
include("DisplacementBoundaryCondition.jl")

# ========================
# Load Boundary Conditions 
# =========================

""" Abstract supertype for all displacement boundary conditions."""
abstract type AbstractLoadBoundaryCondition <: AbstractBoundaryCondition end

" Abstract functor for a `AbstractLoadBoundaryCondition` returns the load at time `t`. "
(lbc::AbstractLoadBoundaryCondition)(t::Real) = values(lbc)(t)

include("GlobalLoadBoundaryCondition.jl")
include("LocalLoadBoundaryCondition.jl")

end # module



