"""
Module defining the boundary conditions implemented.
Two types of boundary conditions are defined Load (Neumann), and Displacements (Dirichlet). 
Overall, each boundary condition consists of a data type with a label, dofs and values into its fields.
"""
module BoundaryConditions

using Reexport: @reexport

@reexport import ..Utils: dofs, label

export AbstractBoundaryCondition, AbstractDisplacementBoundaryCondition, AbstractLoadBoundaryCondition
export DisplacementBoundaryCondition, FixedDofBoundaryCondition, PinnedDisplacementBoundaryCondition
export GlobalLoadBoundaryCondition
# export MᵢLoadBoundaryCondition, MⱼLoadBoundaryCondition, MₖLoadBoundaryCondition
# export FᵢLoadBoundaryCondition, FⱼLoadBoundaryCondition, FₖLoadBoundaryCondition


""" Abstract supertype for all elements.

An `AbstractBoundaryCondition` object facilitates the process of defining:

    - Displacements (Dirichlet) boundary conditions.
    - Load (Neumann) boundary conditions.

**Common methods:**


* [`apply!`](@ref)
* [`label`](@ref)
"""

abstract type AbstractBoundaryCondition end

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

""" Generalized displacement boundary condition struct.
### Fields:
- `dofs`    -- Vectors of symbols the where the boundary condition is subscripted. 
- `values`  -- Values imposed function. 
- `name`    -- Boundary condition label.
"""
Base.@kwdef struct DisplacementBoundaryCondition <: AbstractDisplacementBoundaryCondition
    dofs::Vector{Symbol}
    values::Function
    name::Symbol = :no_labelled_bc
end


""" Fixed displacement boundary condition struct:
This is a particular instance of the struct `DisplacementBoundaryCondition`
    considering null dof value at an specific component of the dof displacements
### Fields:
- `dofs`             -- Vectors of symbols the where the boundary condition is subscribed.
- `fixed_components` -- Vectors of integer indicaating the degree of freedom component fixed.
- `name`             -- Boundary condition label.
"""
Base.@kwdef struct FixedDofBoundaryCondition <: AbstractDisplacementBoundaryCondition
    dofs::Vector{Symbol}
    fixed_components::Vector{Int}
    name::Symbol = :no_labelled_bc
end

FixedDofBoundaryCondition(dofs::Vector{Symbol}, fixed_components::Vector{Int}, name::String="no_labelled_bc") =
    FixedDofBoundaryCondition(dofs, fixed_components, Symbol(name))

fixed_components(bc::FixedDofBoundaryCondition) = bc.fixed_components

# ========================
# Load Boundary Conditions 
# =========================

""" Abstract supertype for all displacement boundary conditions."""

abstract type AbstractLoadBoundaryCondition <: AbstractBoundaryCondition end

""" Load boundary condition imposed in local coordinates of the element.
### Fields:
- `dofs`   -- Degrees of freedom where the boundary condition is imposed. 
- `values` -- Values imposed function. 
- `name`   -- Boundary condition label.
"""
struct LocalLoadBoundaryCondition <: AbstractLoadBoundaryCondition
    dofs::Vector{Symbol}
    values::Function
    name::Symbol
end

""" Load boundary condition imposed in global coordinates of the element.
### Fields:
- `dofs`   -- Degrees of freedom where the boundary condition is imposed. 
- `values` -- Values imposed function. 
- `name`   -- Boundary condition label.
"""
struct GlobalLoadBoundaryCondition <: AbstractLoadBoundaryCondition
    dofs::Vector{Symbol}
    values::Function
    name::Symbol
end

GlobalLoadBoundaryCondition(dofs::Vector{Symbol}, values::Function, name::String="no_labelled_bc") =
    GlobalLoadBoundaryCondition(dofs, values, Symbol(name))

""" Spring boundary condition imposed in global coordinates of the element.
### Fields:
- `dofs`             -- Degrees of freedom where the spring force is imposed. 
- `spring_constants` -- Spring constant `k = f/u` for each dof. 
"""

# ========================
# Spring boundary conditions 
# =========================

struct SpringsBoundaryCondition{D,V} <: AbstractBoundaryCondition
    dofs::D
    spring_constants::V
end


end # module



