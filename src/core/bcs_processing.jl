using ..Utils: row_vector
using ..StructuralModel: AbstractStructure, load_bcs, free_dofs
using ..StructuralAnalyses: AbstractStructuralAnalysis


"Updates the external forces vector with loads boundary conditions at time `t`."
function _update_load_bcs!(s::AbstractStructure, sa::AbstractStructuralAnalysis, t::Number)

    # Extract load boundary conditions
    bcs = boundary_conditions(s)

    for lbc in load_bcs(bcs)

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

        sa.state.Fₑₓₜᵏ[dofs_lbc] = repeat(dofs_values, outer=repeat_mod)
    end
end

"Apply fixed boundary_conditions on the strcture `s`."
function _apply_fixed_bc!(s::AbstractStructure, sa::AbstractStructuralAnalysis)

    # Extract load boundary conditions
    bcs = boundary_conditions(s)
    d_bcs = displacement_bcs(bcs)

    # Filter fixed boundary conditions
    fixed_bcs = filter(bc -> bc isa Union{PinnedDisplacementBoundaryCondition,FixedDisplacementBoundaryCondition}, d_bcs)


    for fbc in fixed_bcs

        # Extract dofs to apply the bc
        fbc_dofs_symbols = dofs(fbc)

        # Extract nodes and elements 
        entities = bcs[fbc]
        dofs_fbc = Dof[]

        for dof_symbol in fbc_dofs_symbols
            dofs_fbc_symbol = row_vector(getindex.(dofs.(entities), dof_symbol))
            push!(dofs_fbc, dofs_fbc_symbol...)
        end

        # Repeat the bc values vector to fill a vector of dofs
        dofs_values = values(fbc)(0.0)
        repeat_mod = Int(length(dofs_fbc) / length(dofs_values))
        sa.state.Uᵏ[dofs_fbc] = repeat(dofs_values, outer=repeat_mod)

        # Delete dofs
        deleteat!(free_dofs(s), findall(x -> x ∈ dofs_fbc, free_dofs(s)))
    end
end

#=
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


=#