using .Elements: dofs
using .BoundaryConditions: AbstractDisplacementBoundaryCondition, FixedDisplacementBoundaryCondition,
    PinnedDisplacementBoundaryCondition, values, load_factor_function
using .StructuralModel: AbstractStructure, sets, load_bcs, disp_bcs, structural_bcs


"Updates the external forces vector with loads boundary conditions at time `t`."
function _apply_load_bc!(s::AbstractStructure, t::Number)

    # Extract load boundary conditions
    s_load_bcs = load_bcs(s)
    bcs_sets = sets(structural_bcs(s))

    for lbc in s_load_bcs
        # Extract load boundary condition info
        lbc_label = label(lbc)
        lbc_dofs_symbols = dofs(lbc)
        lbc_values = values(lbc)

        # Find dofs of the element to apply that load 
        elem_indexes_lbc = bcs_sets[string(lbc_label)]
        elems_lbc = s[elem_indexes_lbc]
        dof_elems_bc = (dofs(elems_lbc))

        # Filter dofs of the element with same symbol 
        for dof_symbol in lbc_dofs_symbols
            lbc_dofs = filter(d -> symbol(d) == dof_symbol, dof_elems_bc)
            s.state.Fₑₓₜᵏ[lbc_dofs] = load_factor_function(lbc)(t) * lbc_values
        end
    end
end


"Updates displacements vector with displacements boundary conditions at time `t`."
function _apply_disp_bc!(s::AbstractStructure, t::Number)

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

    # Extract displacement boundary condition info
    dbc_label = label(dbc)
    dbc_dofs_symbols = dofs(dbc)
    dbc_values = values(dbc)

    # Find dofs of the element to apply that displacement 
    elem_indexes_dbc = bcs_sets[string(dbc_label)]
    elems_dbc = s[elem_indexes_dbc]
    dof_elems_bc = (dofs(elems_dbc))
    # Filter dofs of the element with same symbol 
    for (i, dof_symbol) in enumerate(dbc_dofs_symbols)
        dbc_dofs = filter(d -> symbol(d) == dof_symbol, dof_elems_bc)
        s.state.Uᵏ[dbc_dofs] = dbc_values(t)[i] * ones(length(dbc_dofs))
    end
end

function _fill_displacements_bc!(s::AbstractStructure, dbc::D, bcs_sets::Dict, t::Real) where
{D<:Union{FixedDisplacementBoundaryCondition,PinnedDisplacementBoundaryCondition}}

    # Extract displacement boundary condition info
    dbc_label = label(dbc)
    dbc_dofs_symbols = dofs(dbc)
    dbc_values = values(dbc)

    # Find dofs of the element and fix! the dofs which have the same symbol as the fixed
    elem_indexes_dbc = bcs_sets[string(dbc_label)]
    elems_dbc = s[elem_indexes_dbc]
    dbc_dofs = dofs(elems_dbc)

    for dof_symbol in dbc_dofs_symbols
        dbc_dofs_to_fix = filter(d -> symbol(d) == dof_symbol, dbc_dofs)
        fix!.(dbc_dofs_to_fix)
        s.state.Uᵏ[dbc_dofs_to_fix] = zeros(length(dbc_dofs_to_fix))
    end

end

