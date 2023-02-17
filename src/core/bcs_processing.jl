using ..Utils: row_vector

using ..BoundaryConditions: AbstractLoadBoundaryCondition, AbstractDisplacementBoundaryCondition, FixedDofBoundaryCondition, dofs, fixed_components
using ..StructuralModel: AbstractStructure, load_bcs, free_dofs
using ..StructuralAnalyses: AbstractStructuralAnalysis, external_forces

export _apply!

"Applies a fixed displacement boundary condition to the structure `s` "
function _apply!(s::AbstractStructure, fbc::FixedDofBoundaryCondition)

    bcs = boundary_conditions(s)

    # Extract dofs to apply the bc
    fbc_dofs_symbols = dofs(fbc)

    # Extract nodes and elements 
    entities = bcs[fbc]

    for dof_symbol in fbc_dofs_symbols
        dofs_entities = getindex.(dofs.(entities), dof_symbol)
        for component in fixed_components(fbc)
            dofs_to_delete = getindex.(dofs_entities, component)
            deleteat!(free_dofs(s), findall(x -> x ∈ dofs_to_delete, free_dofs(s)))
        end
    end

end

"Applies a vector of fixed displacement condition to the structure `s` "
_apply!(s::AbstractStructure, l_bcs::Vector{<:FixedDofBoundaryCondition}) = [_apply!(s, lbc) for lbc in l_bcs]


"Applies a fixed displacement boundary condition to the structural analysis `sa` at the current analysis time `t`"
function _apply!(sa::AbstractStructuralAnalysis, lbc::AbstractLoadBoundaryCondition)

    t = current_time(sa)
    bcs = boundary_conditions(structure(sa))

    # Extract dofs to apply the bc
    lbc_dofs_symbols = dofs(lbc)

    # Extract nodes and elements 
    entities = bcs[lbc]
    dofs_lbc = Dof[]

    for dof_symbol in lbc_dofs_symbols
        dofs_lbc_symbol = row_vector(getindex.(dofs.(entities), dof_symbol))
        push!(dofs_lbc, dofs_lbc_symbol...)
    end

    # Repeat the bc values vector to fill a vector of dofs
    dofs_values = values(lbc)(t)

    repeat_mod = Int(length(dofs_lbc) / length(dofs_values))

    external_forces(current_state(sa))[dofs_lbc] = repeat(dofs_values, outer=repeat_mod)

end

"Applies a vector of load boundary conditions to the structure `s` "
_apply!(sa::AbstractStructuralAnalysis, l_bcs::Vector{<:AbstractLoadBoundaryCondition}) = [_apply!(sa, lbc) for lbc in l_bcs]




