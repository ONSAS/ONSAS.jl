"""
Module defining the boundary conditions implemented.
"""
module BoundaryConditions

using Reexport: @reexport

@reexport import ..Utils: dofs, label

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
* [`values`](@ref)
"""

abstract type AbstractBoundaryCondition end

const DEFAULT_LABEL = :bc_label_no_assigned

"Returns the degrees of freedom where the boundary condition is imposed"
dofs(bc::AbstractBoundaryCondition) = bc.dofs

label(bc::AbstractBoundaryCondition) = bc.name

"Returns the values imposed to the respective degrees of freedom"
Base.values(bc::AbstractBoundaryCondition) = bc.values

# ================================
# Displacement Boundary Conditions 
# ================================

""" Abstract supertype for all displacement boundary conditions."""
abstract type AbstractDisplacementBoundaryCondition <: AbstractBoundaryCondition end


""" Generalized displacement boundary condition struct.
### Fields:
- `dofs`    -- Local degrees of freedom where the boundary condition is imposed. 
- `values`  -- Values imposed. 
- `name`    -- Boundary condition label.
"""
Base.@kwdef struct DisplacementBoundaryCondition{D,F} <: AbstractDisplacementBoundaryCondition
    dofs::D
    values::F
    name::Symbol = DEFAULT_LABEL
end

function DisplacementBoundaryCondition(dofs, values::F, label::L) where {F,L}
    DisplacementBoundaryCondition(dofs, values, Symbol(label))
end

""" Fixed displacement boundary condition struct:

This is a particular instance of the struct `DisplacementBoundaryCondition`
    considering null displacements and rotations.
### Fields:
- `bc`    -- Displacement boundary condition constructed with fixed dofs and values. 
"""
struct FixedDisplacementBoundaryCondition <: AbstractDisplacementBoundaryCondition
    bc::DisplacementBoundaryCondition
    function FixedDisplacementBoundaryCondition(label_bc=DEFAULT_LABEL)

        local_dofs_fixed = [1, 2, 3, 4, 5, 6]

        bc = DisplacementBoundaryCondition(
            dofs=local_dofs_fixed,
            values=t -> zeros(length(local_dofs_fixed)),
            name=label_bc
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
    function PinnedDisplacementBoundaryCondition(label_bc=DEFAULT_LABEL)

        local_dofs_fixed = [1, 3, 5]

        bc = DisplacementBoundaryCondition(
            dofs=local_dofs_fixed,
            values=t -> zeros(length(local_dofs_fixed)),
            name=label_bc
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

const DEFAULT_LOAD_FACTOR_FUNC = t -> t

""" Generalized load boundary condition imposed in local coordinates of the element.
### Fields:
- `dofs`        -- Degrees of freedom where the boundary condition is imposed. 
- `values`      -- Values imposed. 
- `load_factor` -- Function to compute at each time step the load factor.
- `name`       -- Boundary condition label.
"""
Base.@kwdef struct LocalLoadBoundaryCondition{D,V,F} <: AbstractLoadBoundaryCondition
    dofs::D
    values::V
    load_factor::F = DEFAULT_LOAD_FACTOR_FUNC
    name::Symbol = DEFAULT_LABEL
end

"Returns the load factor function"
load_factor_function(lbc::LocalLoadBoundaryCondition) = lbc.load_time_factor

""" Load boundary condition imposed in global coordinates of the element.
### Fields:
- `dofs`        -- Degrees of freedom where the boundary condition is imposed. 
- `values`      -- Values imposed. 
- `load_factor` -- Function to compute at each time step the load factor.
- `name`       -- Boundary condition label.
"""
Base.@kwdef struct GlobalLoadBoundaryCondition{D,V,F} <: AbstractLoadBoundaryCondition
    dofs::D
    values::V
    load_factor::F = DEFAULT_LOAD_FACTOR_FUNC
    name::Symbol = DEFAULT_LABEL
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
struct MᵢLoadBoundaryCondition{V,F} <: AbstractLoadBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function MᵢLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {V<:Real,F}
        m_dofs = [2]#θᵢ
        lbc = GlobalLoadBoundaryCondition(m_dofs, [v], load_factor, label)
        return new{V,F}(lbc)
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
struct MⱼLoadBoundaryCondition{V,F} <: AbstractLoadBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function MⱼLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {V<:Real,F}
        m_dofs = [4]
        lbc = GlobalLoadBoundaryCondition(m_dofs, [v], load_factor, label)
        return new{V,F}(lbc)
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
struct MₖLoadBoundaryCondition{V,F} <: AbstractLoadBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function MₖLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {V<:Real,F}
        m_dofs = [6]
        lbc = GlobalLoadBoundaryCondition(m_dofs, [v], load_factor, label)
        return new{V,F}(lbc)
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
struct FᵢLoadBoundaryCondition{V,F} <: AbstractLoadBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function FᵢLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {V<:Real,F}
        m_dofs = [1]
        lbc = GlobalLoadBoundaryCondition(m_dofs, [v], load_factor, label)
        return new{V,F}(lbc)
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
struct FⱼLoadBoundaryCondition{V,F} <: AbstractLoadBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function FⱼLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {V<:Real,F}
        m_dofs = [3]
        lbc = GlobalLoadBoundaryCondition(m_dofs, [v], load_factor, label)
        return new{V,F}(lbc)
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
struct FₖLoadBoundaryCondition{V,F} <: AbstractLoadBoundaryCondition
    bc::GlobalLoadBoundaryCondition
    function FₖLoadBoundaryCondition(
        v::V, load_factor::F=DEFAULT_LOAD_FACTOR_FUNC; label=DEFAULT_LABEL
    ) where {V<:Real,F}
        m_dofs = [5]
        lbc = GlobalLoadBoundaryCondition(m_dofs, [v], load_factor, label)
        return new{V,F}(lbc)
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
struct UserLoadsBoundaryCondition{F} <: AbstractLoadBoundaryCondition
    f::F
end

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



