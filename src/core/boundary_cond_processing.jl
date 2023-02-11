using ..Elements: dofs, index
using ..Utils: dimension, label
using ..BoundaryConditions: AbstractDisplacementBoundaryCondition, DisplacementBoundaryCondition, FixedDisplacementBoundaryCondition,
    PinnedDisplacementBoundaryCondition, values, load_factor_function
using ..StructuralModel: AbstractStructure, mesh, sets, load_bcs, disp_bcs, structural_bcs
using ..StructuralAnalyses: AbstractStructuralAnalysis, AbstractStructuralState, structure


"Updates the external forces vector with loads boundary conditions at time `t`."
function _apply_load_bc!(s::AbstractStructure, sa::AbstractStructuralAnalysis, t::Number)

    # Extract load boundary conditions
    s_load_bcs = load_bcs(s)
    bcs_sets = sets(structural_bcs(s))

    for lbc in s_load_bcs
        # Extract load boundary condition info
        lbc_label = label(lbc)
        lbc_dofs_indexes = dofs(lbc)
        lbc_values = values(lbc)

        # Find dofs of the element to apply that load 
        elem_indexes_lbc = bcs_sets[string(lbc_label)]
        elems_lbc = s[elem_indexes_lbc]
        dof_elems_bc = (dofs(elems_lbc))

        # Filter dofs of the element with same index 
        for lbc_dof_index in lbc_dofs_indexes
            lbc_dofs_indexes = [dof_elems_bc[lbc_dof_index]]
            sa.state.Fₑₓₜᵏ[lbc_dofs_indexes] = load_factor_function(lbc)(t) * lbc_values
        end
    end
end

"Updates displacements vector with displacements boundary conditions at time `t`."
function _apply_disp_bc!(s::AbstractStructure, sa::AbstractStructuralAnalysis, t::Number)

    # Extract disp boundary conditions
    s_displacement_bcs = disp_bcs(s)
    bcs_sets = sets(structural_bcs(s))

    for dbc in s_displacement_bcs
        _fill_displacements_bc!(s, dbc, bcs_sets, t)
    end
end

"Fills the structure displacements vector of the current state."
function _fill_displacements_bc!(
    s::AbstractStructure, dbc::DisplacementBoundaryCondition, bcs_sets::Dict, t::Real)

    error("Please implement non homgenious dirchlet boundary conditions")
end

function _fill_displacements_bc!(s::AbstractStructure, dbc::BC, bcs_sets::Dict, t::Real) where
{BC<:Union{FixedDisplacementBoundaryCondition,PinnedDisplacementBoundaryCondition}}

    D = dimension(mesh(s))
    # Extract displacement boundary condition info
    dbc_label = label(dbc)
    dbc_dof_indexes = dofs(dbc)
    dbc_values = values(dbc)

    # Find dofs of the element and fix! the dofs which have the same symbol as the fixed
    elem_indexes_dbc = bcs_sets[string(dbc_label)]
    elems_dbc = s[elem_indexes_dbc]
    dbc_dofs = dofs(elems_dbc)

    indexes_to_fix = Vector{Int}(undef, 0)

    for dof_index in dbc_dof_indexes
        if D == 2
            dbc_dofs_to_fix_indexes = index.(filter(d -> mod(index(d), 3) == dof_index, dbc_dofs))
        elseif D == 3
            dbc_dofs_to_fix_indexes = index.(filter(d -> mod(index(d), 6) == dof_index, dbc_dofs))
        end
        [push!(indexes_to_fix, d) for d in dbc_dofs_to_fix_indexes]
    end
    deleteat!(free_dof_indexes(s), sort(unique(indexes_to_fix)))
end

