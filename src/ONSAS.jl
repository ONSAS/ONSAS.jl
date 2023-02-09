module ONSAS

# load packages
include("init.jl")

# input modules 
include("interface/input_structs.jl")

# structural modules 
include("interface/model_structs.jl")

# core
include("core/update_structure_state.jl")
include("core/ONSAS_init.jl")
include("core/assembler.jl")

end # module
