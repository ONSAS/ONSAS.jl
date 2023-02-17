using ..Utils: row_vector

using ..BoundaryConditions: AbstractLoadBoundaryCondition, AbstractDisplacementBoundaryCondition, FixedDofBoundaryCondition, dofs, fixed_components
using ..StructuralModel: AbstractStructure, load_bcs, free_dofs
using ..StructuralAnalyses: AbstractStructuralAnalysis, external_forces

export _apply!, fixed_components

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




