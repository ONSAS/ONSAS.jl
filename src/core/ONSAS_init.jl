
using .StructuralAnalyses: StaticAnalysis, structure
using .StructuralSolvers: AbstractSolver

import .StructuralSolvers: init

"Returns load factors vector"
_load_factors_vector(sa::StaticAnalysis, alg::AbstractSolver) =
    LinRange(0, final_time(sa), ceil(Int, final_time(sa) / step_size(alg)))


"Returns the initialized analysis and solution struct. "
function init(sa::StaticAnalysis, alg::AbstractSolver, args...; kwargs...)

    s = structure(sa)
    #TODO : add an optinal load factors vector in kwargs 
    λs = _load_factors_vector(sa, alg)

    # Apply load BC into the global external forces vector 
    _apply_load_bc!(s, first(λs))

    # Build initial external forces vector and apply boundary conditions
    _apply_disp_bc!(s, first(λs))

    # Main.@infiltrate


    return sa
end



#=


function generate_neum(nNodes)

    neumDofs = nodes2dofs(collect(1:nNodes), 6)
    neumDofs[2:2:end] .= 0
    return neumDofs
end


function compute_neum_dofs(Mesh, Boundary_conditions)

    nnodes = size(Mesh.nodal_coords, 1)

    # keep non-zero boundary conds indexes
    boundary_types = sort(unique(Mesh.MGBI_mat[:, 3]))
    boundary_types[1] == 0 && popat!(boundary_types, 1)

    for BCnum in boundary_types

    end


    neum_dofs = generate_neum(nnodes)


    # print("neumdofs", neum_dofs)
    # FG = zeros( 6*numNodes )

    # elementsBC = computeElemsByMEBI( MEBIValsMat[:,3], MEBIVec )

    # for indBC in (1:length(BCsData))
    #     print("indBC", indBC)
    #     for elem in elementsBC[indBC]
    #         print("elem:", elem)
    #         nodesElem = elemNodalConnec[ elem ]
    #         dofs      = nodes2dofs( nodesElem, 6)

    #         if length( BCsData[indBC].NeumannNodalDOFs) >0
    #             dofsNeu = dofs[ BCsData[indBC].NeumannNodalDOFs ]
    #             FG[ dofsNeu ] = FG[ dofsNeu ] + BCsData[indBC].NeumannNodalVals
    #         end

    #         if length( BCsData[indBC].DirichletNodalDOFs) >0
    #             @assert sum( BCsData[indBC].DirichletNodalVals )==0 "only homogeneous dirichlet BC are allowed by now!"
    #             print("dofs:", dofs)
    #             neumDofs[ dofs[ BCsData[indBC].DirichletNodalDOFs ] ] .= 0
    #         end
    #     end
    # end

    # print("neumdofs", neumDofs)

    # neumDofs = sort( unique(neumDofs) )
    # neumDofs[1] == 0 && popfirst!(neumDofs)

    # print("neumdofs", neumDofs)
    # print("listo\n")

    return neum_dofs
end

=#