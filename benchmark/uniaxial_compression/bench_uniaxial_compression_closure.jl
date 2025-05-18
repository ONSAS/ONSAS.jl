using BenchmarkTools, ONSAS, Suppressor

# Utils
include("./../bench_utils.jl")

"Print the time to build the structure, define the analysis and solve the analysis."
function print_times(t_structure, t_problem, t_solve, t_point_eval_handler, t_eval_sol)
    @info "Time to build the structure: $t_structure"
    @info "Time to define the analysis: $t_problem"
    @info "Time to solve the analysis: $t_solve"
    @info "Time to build the point evaluation handler: $t_point_eval_handler"
    @info "Time to evaluate the solution: $t_eval_sol"
end

"Runs the experiment and prints the times."
function run_experiment(build_structure::Function,
        analysis,
        alg::AbstractSolver;
        ms, NSTEPS, N_POINTS_EVAL)
    structure, t_structure, _ = @timed build_structure(; ms)

    problem, t_problem = @timed analysis(structure, NSTEPS = NSTEPS)
    reset!(problem)

    solution, t_solve, _ = @timed solve!(problem, alg)

    ph, t_point_eval_handler,
    _ = @timed PointEvalHandler(structure,
        [Point(rand(3)...)
         for i in 1:N_POINTS_EVAL])

    _, t_eval_sol, _ = @timed displacements(solution, ph)

    return t_structure, t_problem, t_solve, t_point_eval_handler, t_eval_sol
end

# Alg to solve static problems
tols = ConvergenceSettings(; rel_U_tol = 1e-8, rel_res_force_tol = 1e-6, max_iter = 20);
alg = NewtonRaphson(tols);
# Static analysis number of steps to reach the final load factor value
NSTEPS = 8;
# Refinement factor 
ms = 0.5;
# Number of points to evaluate the solution
N_POINTS_EVAL = 100;

# ===============================================================
# Uniaxial compression.
# ===============================================================
example_name = "uniaxial_compression";
example_folder, bench_path = joinpath_example_folder(example_name);
include(bench_path);

# Compilation time 
times_compilation = run_experiment(
    uniaxial_compression_structure, NonLinearStaticAnalysis, alg;
    ms = ms, NSTEPS = NSTEPS, N_POINTS_EVAL = N_POINTS_EVAL);
println("Compiling 🚧:")
print_times(times_compilation...)

# First experiment
times₁ = run_experiment(uniaxial_compression_structure, NonLinearStaticAnalysis, alg;
    ms = ms, NSTEPS = NSTEPS, N_POINTS_EVAL = N_POINTS_EVAL);
println("Experiment 1 🔨:")
print_times(times₁...)

# Second experiment
times₂ = run_experiment(uniaxial_compression_structure, NonLinearStaticAnalysis, alg;
    ms = ms, NSTEPS = NSTEPS, N_POINTS_EVAL = N_POINTS_EVAL);
println("Experiment 2 🔨:")
print_times(times₂...)

delete_files(example_folder, ".msh");
