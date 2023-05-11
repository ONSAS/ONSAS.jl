using BenchmarkTools, ONSAS, Suppressor

include("bench_utils.jl");

# Booleans to include examples into benchmarks 
uniaxial_extension = true;
uniaxial_compression = true;
linear_extension = true;
linear_cylinder_internal_pressure = true;

# Parent BenchmarkGroup to contain our suite.
SUITE = BenchmarkGroup();

# Number of evaluations for each example in the suite.
evals = 2;
# Number of samples for each example in the suite.
samples = 2;

# =======================================
# Static analysis benchmarks
# =======================================
# Alg to solve static problems
tols = ConvergenceSettings(; rel_U_tol=1e-8, rel_res_force_tol=1e-6, max_iter=20);
alg = NewtonRaphson(tols);
# Static analysis number of steps to reach the final load factor value
NSTEPS = 8;

if uniaxial_extension
    # ========================================
    # Uniaxial extension.
    # ========================================
    example_name = "uniaxial_extension"
    SUITE[example_name] = BenchmarkGroup()
    example_folder, bench_path = joinpath_example_folder(example_name)
    include(bench_path)

    # Refinement factors 
    ms_range = [0.5 0.4 0.3 0.2 0.1]

    for ms in ms_range
        local structure
        output = @capture_out begin
            structure = uniaxial_extension_structure(; ms)
        end
        nnodes, nelems = num_nodes(structure), num_elements(structure)
        problem = NonLinearStaticAnalysis(structure; NSTEPS=NSTEPS)
        SUITE[example_name]["solve, ms = $ms, nelems = $nelems, nnodes = $nnodes"] = @benchmarkable solve!($problem,
                                                                                                           $alg) evals = evals samples = samples
    end

    # Remove all .msh files from the example_folder 
    delete_files(example_folder, ".msh")
end

if uniaxial_compression
    # ===============================================================
    # Uniaxial compression.
    # ===============================================================
    example_name = "uniaxial_compression"
    SUITE[example_name] = BenchmarkGroup()
    example_folder, bench_path = joinpath_example_folder(example_name)
    include(bench_path)

    for ms in ms_range
        local structure
        output = @capture_out begin
            structure = uniaxial_compression_structure(; ms)
        end
        # Extract the number of elements as an integer.
        nnodes, nelems = num_nodes(structure), num_elements(structure)
        problem = NonLinearStaticAnalysis(structure; NSTEPS=NSTEPS)
        SUITE[example_name]["solve, ms = $ms, nelems = $nelems, nnodes = $nnodes"] = @benchmarkable solve!($problem,
                                                                                                           $alg) evals = evals samples = samples
    end

    # Remove all .msh files from the example_folder 
    delete_files(example_folder, ".msh")
end

if linear_extension
    # ===============================================================
    # Linear extension.
    # ===============================================================
    example_name = "linear_extension"
    SUITE[example_name] = BenchmarkGroup()
    example_folder, bench_path = joinpath_example_folder(example_name)
    include(bench_path)

    for ms in ms_range
        local structure
        output = @capture_out begin
            structure = linear_extension_structure(; ms)
        end
        # Extract the number of elements as an integer.
        nnodes, nelems = num_nodes(structure), num_elements(structure)
        problem = LinearStaticAnalysis(structure; NSTEPS=NSTEPS)
        SUITE[example_name]["solve, ms = $ms, nelems = $nelems, nnodes = $nnodes"] = @benchmarkable solve!($problem) evals = evals samples = samples
    end

    # Remove all .msh files from the example_folder 
    delete_files(example_folder, ".msh")
end

if linear_cylinder_internal_pressure
    # ===============================================================
    # Linear cylinder with internal pressure.
    # ===============================================================
    example_name = "linear_cylinder_internal_pressure"

    SUITE[example_name] = BenchmarkGroup()
    example_folder, bench_path = joinpath_example_folder(example_name)
    include(bench_path)

    ms_range = [0.5 1.0]
    n_points_each_axis_range = [5, 10]

    for NPOINTS in n_points_each_axis_range
        for ms in ms_range
            local structure
            structure = linear_cylinder_structure(; ms)

            nnodes, nelems = num_nodes(structure), num_elements(structure)
            SUITE[example_name]["structure, ms = $ms, nelems = $nelems, nnodes = $nnodes"] = @benchmarkable linear_cylinder_structure(ms=$ms) evals = evals samples = samples

            problem = LinearStaticAnalysis(structure; NSTEPS)

            solution = solve!(problem)
            reset!(problem)
            SUITE[example_name]["solve!, ms = $ms, nelems = $nelems, nnodes = $nnodes"] = @benchmarkable solve!($problem) evals = evals samples = samples

            ph = point_eval_handler(structure; NPOINTS)
            SUITE[example_name]["point_eval_handler, ms = $ms, nelems = $nelems, nnodes = $nnodes, npoints = $(NPOINTS^3)"] = @benchmarkable point_eval_handler($structure,
                                                                                                                                                                NPOINTS=$NPOINTS) evals = evals samples = samples
            eval_solution = displacements(solution, ph)
        end
    end
    # Remove all .msh files from the example_folder 
    delete_files(example_folder, ".msh")
end
