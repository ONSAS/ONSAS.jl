#################################
# Boundary Conditions interface #
#################################

"""
Module defining the boundary conditions implemented.
"""
module BoundaryConditions


export AbstractBoundaryCondition, label, imposed_dofs, imposed_values


""" Abstract supertype for a boundary condition.
The following methods are provided by the interface:
- `label(bc)`  -- return the boundary condition label 
- `dofs(bc)`   -- return the degrees of freedom where the boundary condition is imposed 
- `values(bc)` -- return the values where the boundary condition is imposed 
"""

abstract type AbstractBoundaryCondition end

const DEFAULT_LABEL = :bc_label_no_assigned

"Returns the degrees of freedom where the boundary condition is imposed"
dofs(bc::ABC) = bc.dofs

"Returns the boundary condition label"
label(bc::ABC) = bc.label

"Returns the values imposed to the respective degrees of freedom"
values(bc::ABC) = bc.values


#####################################
# Displacements Boundary Conditions #
#####################################

""" Displacement boundary condition struct.
### Fields:
- `dofs`    -- Degrees of freedom where the boundary condition is imposed. 
- `vals`    -- Values imposed. 
- `label`   -- Boundary condition label.
"""
Base.@kwdef struct DisplacementBoundaryCondition{D,V} <: AbstractBoundaryCondition
    dofs::D
    values::V
    label::Symbol = DEFAULT_LABEL
end

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

end # module



