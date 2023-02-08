module ONSAS

# load packages
include("init.jl")

# interface
include("interface/input_structs.jl")
include("interface/model_structs.jl")

# core
include("core/update_structure_state.jl")
include("core/ONSAS_init.jl")

end # module
