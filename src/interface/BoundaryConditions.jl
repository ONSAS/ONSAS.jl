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
export GlobalLoadBoundaryCondition, load_factor_function
export MᵢLoadBoundaryCondition, MⱼLoadBoundaryCondition, MₖLoadBoundaryCondition
export FᵢLoadBoundaryCondition, FⱼLoadBoundaryCondition, FₖLoadBoundaryCondition


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

""" Generalized displacement boundary condition struct.
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

dofs(pbc::PinnedDisplacementBoundaryCondition) = dofs(pbc.bc)
Base.values(pbc::PinnedDisplacementBoundaryCondition) = values(pbc.bc)
label(pbc::PinnedDisplacementBoundaryCondition) = label(pbc.bc)

#############################
# Load Boundary Conditions #
#############################

const DEFAULT_LOAD_FACTOR_FUNC = (t) -> 1.0

""" Load boundary condition imposed in local coordinates of the element.
### Fields:
- `dofs`        -- Degrees of freedom where the boundary condition is imposed. 
- `vals`        -- Values imposed. 
- `load_factor` -- Function to compute at each time step the load factor.
- `label`       -- Boundary condition label.
"""
Base.@kwdef struct LocalLoadBoundaryCondition{D,V,F} <: AbstractBoundaryCondition
    dofs::D
    values::V
    load_factor::F = DEFAULT_LOAD_FACTOR_FUNC
    label::Symbol = DEFAULT_LABEL
end

"Returns the load factor function"
load_factor_function(lbc::LocalLoadBoundaryCondition) = lbc.load_time_factor

""" Load boundary condition imposed in global coordinates of the element.
### Fields:
- `dofs`        -- Degrees of freedom where the boundary condition is imposed. 
- `values`      -- Values imposed. 
- `load_factor` -- Function to compute at each time step the load factor.
- `label`       -- Boundary condition label.
"""
Base.@kwdef struct GlobalLoadBoundaryCondition{D,V,F} <: AbstractBoundaryCondition
    dofs::D
    values::V
    load_factor::F = DEFAULT_LOAD_FACTOR_FUNC
    label::Symbol = DEFAULT_LABEL
end

load_factor_function(gbc::GlobalLoadBoundaryCondition) = gbc.load_factor

""" Mᵢ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive moment along the `i` axis.
### Fields:
- `bc`   -- Mᵢ Load boundary condition. 
"""
struct MᵢLoadBoundaryCondition{F,V} <: AbstractBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function MᵢLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {F,V}
        m_dofs = [:Mᵢ]
        lbc = GlobalLoadBoundaryCondition(m_dofs, v, load_factor, label)
        return new{F,V}(lbc)
    end
end

dofs(Mᵢbc::MᵢLoadBoundaryCondition) = dofs(Mᵢbc.bc)
Base.values(Mᵢbc::MᵢLoadBoundaryCondition) = values(Mᵢbc.bc)
label(Mᵢbc::MᵢLoadBoundaryCondition) = label(Mᵢbc.bc)
load_factor_function(Mᵢbc::MᵢLoadBoundaryCondition) = load_factor_function(Mᵢbc.bc)

""" Mⱼ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive moment along the `j` axis.
### Fields:
- `bc`   -- Mⱼ Load boundary condition. 
"""
struct MⱼLoadBoundaryCondition{F,V} <: AbstractBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function MⱼLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {F,V}
        m_dofs = [:Mⱼ]
        lbc = GlobalLoadBoundaryCondition(m_dofs, v, load_factor, label)
        return new{F,V}(lbc)
    end
end

dofs(Mⱼbc::MⱼLoadBoundaryCondition) = dofs(Mⱼbc.bc)
Base.values(Mⱼbc::MⱼLoadBoundaryCondition) = values(Mⱼbc.bc)
label(Mⱼbc::MⱼLoadBoundaryCondition) = label(Mⱼbc.bc)
load_factor_function(Mⱼbc::MⱼLoadBoundaryCondition) = load_factor_function(Mⱼbc.bc)

""" Mₖ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive moment along the `j` axis.
### Fields:
- `bc`   -- Mₖ Load boundary condition. 
"""
struct MₖLoadBoundaryCondition{F,V} <: AbstractBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function MₖLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {F,V}
        m_dofs = [:Mₖ]
        lbc = GlobalLoadBoundaryCondition(m_dofs, v, load_factor, label)
        return new{F,V}(lbc)
    end
end

dofs(Mₖbc::MₖLoadBoundaryCondition) = dofs(Mₖbc.bc)
Base.values(Mₖbc::MₖLoadBoundaryCondition) = values(Mₖbc.bc)
label(Mₖbc::MₖLoadBoundaryCondition) = label(Mₖbc.bc)
load_factor_function(Mₖbc::MₖLoadBoundaryCondition) = load_factor_function(Mₖbc.bc)

""" Fᵢ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive force along the `i` axis.
### Fields:
- `bc`  -- Fᵢ Load boundary condition. 
"""
struct FᵢLoadBoundaryCondition{F,V} <: AbstractBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function FᵢLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {F,V}
        m_dofs = [:Fᵢ]
        lbc = GlobalLoadBoundaryCondition(m_dofs, v, load_factor, label)
        return new{F,V}(lbc)
    end
end

dofs(Fᵢbc::FᵢLoadBoundaryCondition) = dofs(Fᵢbc.bc)
Base.values(Fᵢbc::FᵢLoadBoundaryCondition) = values(Fᵢbc.bc)
label(Fᵢbc::FᵢLoadBoundaryCondition) = label(Fᵢbc.bc)
load_factor_function(Fᵢbc::FᵢLoadBoundaryCondition) = load_factor_function(Fᵢbc.bc)


""" Fⱼ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive force along the `i` axis.
### Fields:
- `bc`  -- Fⱼ Load boundary condition. 
"""
struct FⱼLoadBoundaryCondition{F,V} <: AbstractBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function FⱼLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {F,V}
        m_dofs = [:Fⱼ]
        lbc = GlobalLoadBoundaryCondition(m_dofs, v, load_factor, label)
        return new{F,V}(lbc)
    end
end

dofs(Fⱼbc::FⱼLoadBoundaryCondition) = dofs(Fⱼbc.bc)
Base.values(Fⱼbc::FⱼLoadBoundaryCondition) = values(Fⱼbc.bc)
label(Fⱼbc::FⱼLoadBoundaryCondition) = label(Fⱼbc.bc)
load_factor_function(Fⱼbc::FⱼLoadBoundaryCondition) = load_factor_function(Fⱼbc.bc)


""" Fₖ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive force along the `i` axis.
### Fields:
- `bc`  -- Fₖ Load boundary condition. 
"""
struct FₖLoadBoundaryCondition{F,V} <: AbstractBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function FₖLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {F,V}
        m_dofs = [:Fₖ]
        lbc = GlobalLoadBoundaryCondition(m_dofs, v, load_factor, label)
        return new{F,V}(lbc)
    end
end

dofs(Fₖbc::FₖLoadBoundaryCondition) = dofs(Fₖbc.bc)
Base.values(Fₖbc::FₖLoadBoundaryCondition) = values(Fₖbc.bc)
label(Fₖbc::FₖLoadBoundaryCondition) = label(Fₖbc.bc)
load_factor_function(Fₖbc::FₖLoadBoundaryCondition) = load_factor_function(Fₖbc.bc)

""" User Load boundary condition imposed in global coordinates to the entire structure.
### Fields:
- `f`    -- User load function (the output of this function must be a vector with length equals to the number of dofs). 
"""
struct UserLoadsBoundaryCondition <: AbstractBoundaryCondition
    f::Function
end

""" Spring boundary condition imposed in global coordinates of the element.
### Fields:
- `dofs`             -- Degrees of freedom where the spring force is imposed. 
- `spring_constants` -- Spring constant `k = f/u` for each dof. 
"""

struct SpringsBoundaryCondition{D,V} <: AbstractBoundaryCondition
    dofs::D
    spring_constants::V
end


end # module



