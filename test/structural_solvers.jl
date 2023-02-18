##########################
# Structural model tests #
##########################
using Test: @testset, @test
using ONSAS.StructuralSolvers
using LinearAlgebra: norm
using SparseArrays: SparseMatrixCSC, sparse


@testset "ONSAS.StructuralSolvers.Criteria" begin

    @test ResidualForceCriterion <: AbstractConvergenceCriterion
    @test ΔUCriterion <: AbstractConvergenceCriterion
    @test ΔU_and_ResidualForce_Criteria <: AbstractConvergenceCriterion
    @test MaxIterCriterion <: AbstractConvergenceCriterion
    @test NotConvergedYet <: AbstractConvergenceCriterion

end


# Convergence settings
tol_f = 1e-3
tol_u = 1e-5
max_iter = 100
tols = ConvergenceSettings(tol_u, tol_f, max_iter)

@testset "ONSAS.StructuralSolvers.ResidualsIterationStep" begin


    @test residual_forces_tol(tols) == tol_f
    @test displacement_tol(tols) == tol_u
    @test max_iter_tol(tols) == max_iter

    # Residuals
    residuals_current_step = ResidualsIterationStep()

    @test iter(residuals_current_step) == 0
    @test criterion(residuals_current_step) isa NotConvergedYet

    _reset!(residuals_current_step)
    @test iter(residuals_current_step) == 0
    @test criterion(residuals_current_step) isa NotConvergedYet
    @test all(map(x -> x ≥ 1e3, displacement_tol(residuals_current_step)))
    @test all(map(x -> x ≥ 1e3, residual_forces_tol(residuals_current_step)))
    _update!(residuals_current_step, ΔUCriterion())
    @test !isconverged!(residuals_current_step, tols)

    # Update residuals
    ΔU = [1e-10, 1e-10]
    U = [1e-1, 1e-4]
    Δr = [1e-10, 1e-10]
    fₑₓₜ = [1e3, 1e3]

    _update!(residuals_current_step, ΔU, U, Δr, fₑₓₜ)

    @test iter(residuals_current_step) == 1
    ΔU_rel = norm(ΔU) / norm(U)
    Δr_rel = norm(Δr) / norm(fₑₓₜ)
    @test displacement_tol(residuals_current_step)[1] == ΔU_rel
    @test displacement_tol(residuals_current_step)[2] == norm(ΔU)
    @test residual_forces_tol(residuals_current_step)[1] == Δr_rel
    @test residual_forces_tol(residuals_current_step)[2] == norm(Δr)
    @test isconverged!(residuals_current_step, tols)

    # Reset again
    _reset!(residuals_current_step)
    @test iter(residuals_current_step) == 0
    @test criterion(residuals_current_step) isa NotConvergedYet
    @test !isconverged!(residuals_current_step, tols)

end

@testset "ONSAS.StructuralSolvers.NewtonRaphson" begin

    nr = NewtonRaphson(tols)
    tolerances(nr) == tols

end

@testset "ONSAS.StructuralSolvers.Assembler" begin

    Ke = [
        1.0 2.0
        3.0 4.0
    ]

    K_glob = [
        1.0 2.0 0.0
        3.0 5.0 2.0
        0.0 3.0 4.0
    ]

    dofs_element_1 = [Dof(1), Dof(2)]
    dofs_element_2 = [Dof(2), Dof(3)]
    dofs_elements = [dofs_element_1, dofs_element_2]

    I = Int[]
    J = Int[]
    V = Float64[]

    a = Assembler(length(Ke))
    for n_e in 1:2
        _assemble!(a, dofs_elements[n_e], Ke)
    end

    @test a.I == [1, 2, 1, 2, 2, 3, 2, 3]
    @test a.J == [1, 1, 2, 2, 2, 2, 3, 3]
    @test a.V == [
        # Element 1
        Ke[1, 1],
        Ke[2, 1],
        Ke[1, 2],
        Ke[2, 2],
        # Element 2
        Ke[1, 1],
        Ke[2, 1],
        Ke[1, 2],
        Ke[2, 2],
    ]

    K_glob_assembler = end_assemble(a)

    K_to_fill_assembler = SparseMatrixCSC(zeros(3, 3))

    end_assemble!(K_to_fill_assembler, a)

    @test all([K_glob_assembler[ind] == val for (ind, val) in enumerate(K_glob_assembler)])
    @test all([K_glob_assembler[ind] == val for (ind, val) in enumerate(K_to_fill_assembler)])

    _reset!(a)
    @test isempty(a.I)
    @test isempty(a.J)
    @test isempty(a.V)

end
