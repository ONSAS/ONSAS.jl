using BenchmarkTools, ONSAS

# Parent BenchmarkGroup to contain our suite.
SUITE = BenchmarkGroup()

# ===============================================================
# Uniaxial extension.
# ===============================================================
include("uniaxial_extension/uniaxial_extension.jl")

SUITE["Uniaxial extension"] = BenchmarkGroup()

tols = ConvergenceSettings(rel_U_tol=1e-8, rel_res_force_tol=1e-6, max_iter=20)
alg = NewtonRaphson(tols)
structure = uniaxial_extension_structure(ms=0.5)
problem = NonLinearStaticAnalysis(structure, NSTEPS=8)
SUITE["Uniaxial extension"]["256 elements"] = @benchmarkable solve($problem, $alg)
