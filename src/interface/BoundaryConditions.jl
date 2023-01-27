#################################
# Boundary Conditions interface #
#################################

"""
Module defining the boundary conditions implemented.
"""
module BoundaryConditions

using Reexport: @reexport

@reexport import ..Utils: label

export AbstractBoundaryCondition, dofs, label
export DisplacementBoundaryCondition, FixedDisplacementBoundaryCondition, PinnedDisplacementBoundaryCondition


""" Abstract supertype for a boundary condition.
The following methods are provided by the interface:
- `label(bc)`  -- return the boundary condition label 
- `dofs(bc)`   -- return the degrees of freedom where the boundary condition is imposed 
- `values(bc)` -- return the values where the boundary condition is imposed 
"""

abstract type AbstractBoundaryCondition end

const DEFAULT_LABEL = :bc_label_no_assigned

"Returns the degrees of freedom where the boundary condition is imposed"
dofs(bc::AbstractBoundaryCondition) = bc.dofs

"Returns the boundary condition label"
label(bc::AbstractBoundaryCondition) = bc.label

"Returns the values imposed to the respective degrees of freedom"
Base.values(bc::AbstractBoundaryCondition) = bc.values


#####################################
# Displacements Boundary Conditions #
#####################################

""" General displacement boundary condition struct.
### Fields:
- `dofs`    -- Degrees of freedom where the boundary condition is imposed. 
- `values`  -- Values imposed. 
- `label`   -- Boundary condition label.
"""
Base.@kwdef struct DisplacementBoundaryCondition{D,V} <: AbstractBoundaryCondition
    dofs::D
    values::V
    label::Symbol = DEFAULT_LABEL
end

""" Fixed displacement boundary condition struct:

This is a particular instance of the struct `DisplacementBoundaryCondition`
    considering null displacements and rotations.
### Fields:
- `bc`    -- Displacement boundary condition constructed with fixed dofs and values. 
"""
struct FixedDisplacementBoundaryCondition{D,V} <: AbstractBoundaryCondition
    bc::DisplacementBoundaryCondition{D,V}
    function FixedDisplacementBoundaryCondition(label_bc=DEFAULT_LABEL)

        dofs_fixed = [:uᵢ, :θᵢ, :uⱼ, :θⱼ, :uₖ, :θₖ]

        bc = DisplacementBoundaryCondition(
            dofs=dofs_fixed,
            values=zeros(length(dofs_fixed)),
            label=label_bc
        )

        return new{typeof(dofs(bc)),typeof(values(bc))}(bc)
    end
end

# Composed accessors 
#TODO:fixme
# for m in methodswith(AbstractBoundaryCondition)
#     print(m)
#     @eval $m(fbc::FixedDisplacementBoundaryCondition) = $m(fbc.bc)
# end
dofs(fbc::FixedDisplacementBoundaryCondition) = dofs(fbc.bc)
Base.values(fbc::FixedDisplacementBoundaryCondition) = values(fbc.bc)
label(fbc::FixedDisplacementBoundaryCondition) = label(fbc.bc)

""" Pinned displacement boundary condition struct:

This is a particular instance of the struct `DisplacementBoundaryCondition`
    considering null displacements.
### Fields:
- `bc`    -- Displacement boundary condition constructed with pinned dofs and values. 
"""
struct PinnedDisplacementBoundaryCondition{D,V} <: AbstractBoundaryCondition
    bc::DisplacementBoundaryCondition{D,V}
    function PinnedDisplacementBoundaryCondition(label_bc=DEFAULT_LABEL)

        dofs_fixed = [:uᵢ, :uⱼ, :uₖ]

        bc = DisplacementBoundaryCondition(
            dofs=dofs_fixed,
            values=zeros(length(dofs_fixed)),
            label=label_bc
        )

        return new{typeof(dofs(bc)),typeof(values(bc))}(bc)
    end
end

# Composed accessors 
#TODO:fixme
# for m in methodswith(AbstractBoundaryCondition)
#     print(m)
#     @eval $m(fbc::FixedDisplacementBoundaryCondition) = $m(fbc.bc)
# end
dofs(pbc::PinnedDisplacementBoundaryCondition) = dofs(pbc.bc)
Base.values(pbc::PinnedDisplacementBoundaryCondition) = values(pbc.bc)
label(pbc::PinnedDisplacementBoundaryCondition) = label(pbc.bc)

#=

#############################
# Load Boundary Conditions #
#############################

""" Load boundary condition imposed in local coordinates of the element.
### Fields:
- `dofs`    -- Degrees of freedom where the boundary condition is imposed. 
- `vals`    -- Values imposed. 
- `label`   -- Boundary condition label.
"""
Base.@kwdef struct LocalLoadBoundaryCondition{D,V} <: AbstractBoundaryCondition
    dofs::D
    values::V
    load_time_factor::Function = (t) -> 1.0
    label::Symbol = DEFAULT_LABEL
end

""" Load boundary condition imposed in global coordinates of the element.
### Fields:
- `dofs`             -- Degrees of freedom where the boundary condition is imposed. 
- `vals`             -- Values imposed. 
- `label`            -- Boundary condition label.
- `load_time_factor` -- Function to compute the time factor of the load.
"""
Base.@kwdef struct GlobalLoadBoundaryCondition{D,V} <: AbstractBoundaryCondition
    dofs::D
    values::V
    label::Symbol = DEFAULT_LABEL
    load_time_factor::Function = (t) -> 1.0
end

""" User Load boundary condition imposed in global coordinates to entire structure.
### Fields:
- `f`    -- User load function (the output of this function must be a vector with length equals to the number of dofs). 
"""
struct UserLoadsBoundaryCondition <: AbstractBoundaryCondition
    f::Function
end


struct SpringsBoundaryCondition <: AbstractDofs
    # to do
end

=#
end # module



