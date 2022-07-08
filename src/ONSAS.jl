module FEMAssembler

    # load packages
    include("init.jl")

    # initialization of structures
    include("interface.jl")

    # elements
    include("elements/linear_truss.jl")
    include("elements/linear_tetrahedron.jl")

    include("mesh/dofs_computations.jl")
    #include("mesh/mshRead.jl")
    include("assembler.jl")

    # exports
    include("exports.jl")

end # module
