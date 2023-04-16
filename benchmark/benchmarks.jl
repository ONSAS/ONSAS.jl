using BenchmarkTools, ONSAS, Suppressor

# Parent BenchmarkGroup to contain our suite.
SUITE = BenchmarkGroup()

# ===============================================================
# Uniaxial extension.
# ===============================================================
include("uniaxial_extension/uniaxial_extension.jl")

SUITE["Uniaxial extension"] = BenchmarkGroup()

tols = ConvergenceSettings(rel_U_tol=1e-8, rel_res_force_tol=1e-6, max_iter=20)
alg = NewtonRaphson(tols)

for ms in [0.5, 0.4, 0.3, 0.2]
    local structure
    output = @capture_out begin
        structure = uniaxial_extension_structure(; ms)
    end
    # Extract the number of elements as an integer.
    nelems = parse(Int, match(r"(\d+) nodes (\d+) elements", output).captures[2])
    problem = NonLinearStaticAnalysis(structure, NSTEPS=8)
    SUITE["Uniaxial extension"]["solve, ms = $ms, nelems = $nelems"] = @benchmarkable solve($problem, $alg) evals = 3 samples = 2
end
