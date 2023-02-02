#################################
# Boundary Conditions interface #
#################################

"""
Module defining the boundary conditions implemented.
"""
module BoundaryConditions

using Reexport: @reexport

@reexport import ..Utils: ScalarWrapper, dofs, label, set_label!


export AbstractBoundaryCondition, AbstractDisplacementBoundaryCondition, AbstractLoadBoundaryCondition
export DisplacementBoundaryCondition, FixedDisplacementBoundaryCondition, PinnedDisplacementBoundaryCondition
export GlobalLoadBoundaryCondition, load_factor_function
export MᵢLoadBoundaryCondition, MⱼLoadBoundaryCondition, MₖLoadBoundaryCondition
export FᵢLoadBoundaryCondition, FⱼLoadBoundaryCondition, FₖLoadBoundaryCondition


""" Abstract supertype for all elements.

An `AbstractBoundaryCondition` object facilitates the process of defining:

    - Load (Neumann) boundary conditions.
    - Displacements (Dirichlet) boundary conditions.

**Common methods:**

* [`label`](@ref)
* [`set_label!`](@ref)
* [`values`](@ref)
"""

abstract type AbstractBoundaryCondition end

const DEFAULT_LABEL = :bc_label_no_assigned

"Returns the degrees of freedom where the boundary condition is imposed"
dofs(bc::AbstractBoundaryCondition) = bc.dofs

label(bc::AbstractBoundaryCondition) = bc.label[]
set_label!(bc::AbstractBoundaryCondition, label) = bc.label[] = Symbol(label)

"Returns the values imposed to the respective degrees of freedom"
Base.values(bc::AbstractBoundaryCondition) = bc.values


#####################################
# Displacements Boundary Conditions #
#####################################

abstract type AbstractDisplacementBoundaryCondition <: AbstractBoundaryCondition end


""" Generalized displacement boundary condition struct.
### Fields:
- `dofs`    -- Degrees of freedom where the boundary condition is imposed. 
- `values`  -- Values imposed. 
- `label`   -- Boundary condition label.
"""
Base.@kwdef struct DisplacementBoundaryCondition{D,V} <: AbstractDisplacementBoundaryCondition
    dofs::D
    values::V
    label::ScalarWrapper{Symbol} = ScalarWrapper(DEFAULT_LABEL)
end

function DisplacementBoundaryCondition(dofs::D, values::V, label::Symbol) where {D,V}
    DisplacementBoundaryCondition(dofs, values, ScalarWrapper(label))
end
function DisplacementBoundaryCondition(dofs::D, values::V, label::String) where {D,V}
    DisplacementBoundaryCondition(dofs, values, ScalarWrapper(Symbol(label)))
end

""" Fixed displacement boundary condition struct:

This is a particular instance of the struct `DisplacementBoundaryCondition`
    considering null displacements and rotations.
### Fields:
- `bc`    -- Displacement boundary condition constructed with fixed dofs and values. 
"""
struct FixedDisplacementBoundaryCondition{D,V} <: AbstractDisplacementBoundaryCondition
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
set_label!(fbc::FixedDisplacementBoundaryCondition, label) = set_label!(fbc.bc, label)

""" Pinned displacement boundary condition struct:

This is a particular instance of the struct `DisplacementBoundaryCondition`
    considering null displacements.
### Fields:
- `bc`    -- Displacement boundary condition constructed with pinned dofs and values. 
"""
struct PinnedDisplacementBoundaryCondition{D,V} <: AbstractDisplacementBoundaryCondition
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
set_label!(pbc::PinnedDisplacementBoundaryCondition, label) = set_label!(pbc.bc, label)

#############################
# Load Boundary Conditions #
#############################

abstract type AbstractLoadBoundaryCondition <: AbstractBoundaryCondition end


const DEFAULT_LOAD_FACTOR_FUNC = (t) -> 1.0

""" Generalized load boundary condition imposed in local coordinates of the element.
### Fields:
- `dofs`        -- Degrees of freedom where the boundary condition is imposed. 
- `values`      -- Values imposed. 
- `load_factor` -- Function to compute at each time step the load factor.
- `label`       -- Boundary condition label.
"""
Base.@kwdef struct LocalLoadBoundaryCondition{D,V,F} <: AbstractLoadBoundaryCondition
    dofs::D
    values::V
    load_factor::F = DEFAULT_LOAD_FACTOR_FUNC
    label::ScalarWrapper{Symbol} = ScalarWrapper(DEFAULT_LABEL)
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
Base.@kwdef struct GlobalLoadBoundaryCondition{D,V,F} <: AbstractLoadBoundaryCondition
    dofs::D
    values::V
    load_factor::F = DEFAULT_LOAD_FACTOR_FUNC
    label::ScalarWrapper{Symbol} = ScalarWrapper(DEFAULT_LABEL)
end

function GlobalLoadBoundaryCondition(
    dofs, values, load_factor_function, label::L) where {L<:Union{String,Symbol}}
    GlobalLoadBoundaryCondition(dofs, values, load_factor_function, ScalarWrapper(Symbol(label)))
end


load_factor_function(gbc::GlobalLoadBoundaryCondition) = gbc.load_factor

""" Mᵢ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive moment along the `i` axis.
### Fields:
- `bc`   -- Mᵢ Load boundary condition. 
"""
struct MᵢLoadBoundaryCondition{F,V} <: AbstractLoadBoundaryCondition
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
set_label!(Mᵢbc::MᵢLoadBoundaryCondition, label) = set_label!(Mᵢbc.bc, label)
load_factor_function(Mᵢbc::MᵢLoadBoundaryCondition) = load_factor_function(Mᵢbc.bc)

""" Mⱼ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive moment along the `j` axis.
### Fields:
- `bc`   -- Mⱼ Load boundary condition. 
"""
struct MⱼLoadBoundaryCondition{F,V} <: AbstractLoadBoundaryCondition
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
set_label!(Mⱼbc::MⱼLoadBoundaryCondition, label) = set_label!(Mⱼbc.bc, label)
load_factor_function(Mⱼbc::MⱼLoadBoundaryCondition) = load_factor_function(Mⱼbc.bc)

""" Mₖ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive moment along the `j` axis.
### Fields:
- `bc`   -- Mₖ Load boundary condition. 
"""
struct MₖLoadBoundaryCondition{F,V} <: AbstractLoadBoundaryCondition
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
set_label!(Mₖbc::MₖLoadBoundaryCondition, label) = set_label!(Mₖbc.bc, label)
load_factor_function(Mₖbc::MₖLoadBoundaryCondition) = load_factor_function(Mₖbc.bc)

""" Fᵢ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive force along the `i` axis.
### Fields:
- `bc`  -- Fᵢ Load boundary condition. 
"""
struct FᵢLoadBoundaryCondition{F,V} <: AbstractLoadBoundaryCondition
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
set_label!(Fᵢbc::FᵢLoadBoundaryCondition, label) = set_label!(Fᵢbc.bc, label)
load_factor_function(Fᵢbc::FᵢLoadBoundaryCondition) = load_factor_function(Fᵢbc.bc)


""" Fⱼ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive force along the `i` axis.
### Fields:
- `bc`  -- Fⱼ Load boundary condition. 
"""
struct FⱼLoadBoundaryCondition{F,V} <: AbstractLoadBoundaryCondition
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
set_label!(Fⱼbc::FⱼLoadBoundaryCondition, label) = set_label!(Fⱼbc.bc, label)
load_factor_function(Fⱼbc::FⱼLoadBoundaryCondition) = load_factor_function(Fⱼbc.bc)


""" Fₖ load boundary condition struct:

This is a particular instance of the struct `GlobalLoadBoundaryCondition`
    considering a positive force along the `i` axis.
### Fields:
- `bc`  -- Fₖ Load boundary condition. 
"""
struct FₖLoadBoundaryCondition{F,V} <: AbstractLoadBoundaryCondition
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
set_label!(Fₖbc::FₖLoadBoundaryCondition, label) = set_label!(Fₖbc.bc, label)
load_factor_function(Fₖbc::FₖLoadBoundaryCondition) = load_factor_function(Fₖbc.bc)

""" User Load boundary condition imposed in global coordinates to the entire structure.
### Fields:
- `f`    -- User load function (the output of this function must be a vector with length equals to the number of dofs). 
"""
struct UserLoadsBoundaryCondition <: AbstractLoadBoundaryCondition
    f::Function
end

""" Spring boundary condition imposed in global coordinates of the element.
### Fields:
- `dofs`             -- Degrees of freedom where the spring force is imposed. 
- `spring_constants` -- Spring constant `k = f/u` for each dof. 
"""

###############################
# Springs Boundary Conditions #
###############################

struct SpringsBoundaryCondition{D,V} <: AbstractBoundaryCondition
    dofs::D
    spring_constants::V
end


end # module



