module ONSAS

    # load packages
    include("init.jl")

    # initialization of structures
    include("interface.jl")

    # elements
    include("elements/linear_truss.jl")
    include("elements/linear_tetrahedron.jl")

    include("mesh/dofs_computations.jl")

    include("core/ONSAS_init.jl")
    include("core/assembler.jl")

    # exports
    include("exports.jl")

end # module
