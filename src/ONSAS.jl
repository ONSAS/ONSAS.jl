module ONSAS

# load packages
include("init.jl")

# interface
include("interface/input_structs.jl")
include("interface/model_structs.jl")

# elements
include("core/time_step_iteration.jl")
include("core/ONSAS_init.jl")
include("core/ONSAS_solve.jl")
include("core/assembler.jl")
include("core/compute_fext.jl")
include("core/boundary_cond_processing.jl")

# elements
include("elements/linear_truss.jl")
include("elements/linear_tetrahedron.jl")
include("elements/elem2NodalLoads.jl")

# mesh
include("mesh/dofs_computations.jl")


# exports
include("exports.jl")

end # module
