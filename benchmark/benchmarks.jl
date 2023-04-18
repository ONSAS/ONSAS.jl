using BenchmarkTools, ONSAS, Suppressor

include("bench_utils.jl")

# Parent BenchmarkGroup to contain our suite.
SUITE = BenchmarkGroup()

# Number of evaluations for each example in the suite.
evals = 3
# Number of samples for each example in the suite.
samples = 2

# =======================================
# Static analysis benchmarks
# =======================================
# Alg to solve static problems
tols = ConvergenceSettings(rel_U_tol=1e-8, rel_res_force_tol=1e-6, max_iter=20)
alg = NewtonRaphson(tols)
# Static analysis number of steps to reach the final load factor value
NSTEPS = 8

# ========================================
# Uniaxial extension.
# ========================================
example_name = "uniaxial_extension"
SUITE[example_name] = BenchmarkGroup()
example_folder, bench_path = joinpath_example_folder(example_name)
include(bench_path)

# Refinement factors 
ms_range = [0.5 0.4 0.3 0.2]

for ms in ms_range
    local structure
    output = @capture_out begin
        structure = uniaxial_extension_structure(; ms)
    end
    nnodes, nelems = capture_print(output)
    problem = NonLinearStaticAnalysis(structure, NSTEPS=NSTEPS)
    SUITE[example_name]["solve, ms = $ms, nelems = $nelems, nnodes = $nnodes"] =
        @benchmarkable solve!($problem, $alg) evals = evals samples = samples
end

# Remove all .msh files from the example_folder 
delete_files(example_folder, ".msh")

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
    nnodes, nelems = capture_print(output)
    problem = NonLinearStaticAnalysis(structure, NSTEPS=NSTEPS)
    SUITE[example_name]["solve, ms = $ms, nelems = $nelems, nnodes = $nnodes"] =
        @benchmarkable solve!($problem, $alg) evals = evals samples = samples

end

# Remove all .msh files from the example_folder 
delete_files(example_folder, ".msh")

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
    nnodes, nelems = capture_print(output)
    problem = LinearStaticAnalysis(structure, NSTEPS=NSTEPS)
    SUITE[example_name]["solve, ms = $ms, nelems = $nelems, nnodes = $nnodes"] =
        @benchmarkable solve!($problem, $alg) evals = evals samples = samples
end

# Remove all .msh files from the example_folder 
delete_files(example_folder, ".msh")
