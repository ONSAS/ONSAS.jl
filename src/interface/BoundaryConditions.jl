"""
Module defining the boundary conditions implemented.
"""
module BoundaryConditions

using Reexport: @reexport

@reexport import ..Utils: dofs, label

export AbstractBoundaryCondition, AbstractDisplacementBoundaryCondition, AbstractLoadBoundaryCondition
export DisplacementBoundaryCondition, FixedDisplacementBoundaryCondition, PinnedDisplacementBoundaryCondition
export GlobalLoadBoundaryCondition
# export MᵢLoadBoundaryCondition, MⱼLoadBoundaryCondition, MₖLoadBoundaryCondition
# export FᵢLoadBoundaryCondition, FⱼLoadBoundaryCondition, FₖLoadBoundaryCondition


""" Abstract supertype for all elements.

An `AbstractBoundaryCondition` object facilitates the process of defining:

    - Displacements (Dirichlet) boundary conditions.
    - Load (Neumann) boundary conditions.

**Common methods:**

* [`label`](@ref)
* [`values`](@ref)
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
- `dofs`    -- Vectors of symbols where the where the boundary condition is subscripted. 
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
    considering null displacements and rotations.
### Fields:
- `bc`    -- Displacement boundary condition constructed with fixed dofs and values. 
"""
struct FixedDisplacementBoundaryCondition <: AbstractDisplacementBoundaryCondition
    bc::DisplacementBoundaryCondition
    function FixedDisplacementBoundaryCondition(dof_dim::Int=3, label_bc=:no_labelled_bc)

        local_dofs_fixed = [:u, :θ]

        bc = DisplacementBoundaryCondition(
            dofs=local_dofs_fixed,
            values=t -> zeros(dof_dim * length(local_dofs_fixed)),
            name=Symbol(label_bc)
        )

        return new(bc)
    end
end

dofs(fbc::FixedDisplacementBoundaryCondition) = dofs(fbc.bc)
Base.values(fbc::FixedDisplacementBoundaryCondition) = values(fbc.bc)
label(fbc::FixedDisplacementBoundaryCondition) = label(fbc.bc)

""" Pinned displacement boundary condition struct:
This is a particular instance of the struct `DisplacementBoundaryCondition`
    considering null displacements.
### Fields:
- `bc`    -- Displacement boundary condition constructed with pinned dofs and values. 
"""
struct PinnedDisplacementBoundaryCondition <: AbstractDisplacementBoundaryCondition
    bc::DisplacementBoundaryCondition
    function PinnedDisplacementBoundaryCondition(dof_dim::Int=3, label_bc=:no_labelled_bc)

        local_dofs_fixed = [:u]

        bc = DisplacementBoundaryCondition(
            dofs=local_dofs_fixed,
            values=t -> zeros(length(dof_dim * length(local_dofs_fixed))),
            name=Symbol(label_bc)
        )

        return new(bc)
    end
end

dofs(pbc::PinnedDisplacementBoundaryCondition) = dofs(pbc.bc)
Base.values(pbc::PinnedDisplacementBoundaryCondition) = values(pbc.bc)
label(pbc::PinnedDisplacementBoundaryCondition) = label(pbc.bc)


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



