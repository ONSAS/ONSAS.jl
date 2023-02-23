using ..Meshes: AbstractMesh
using ..BoundaryConditions: FixedDofBoundaryCondition
using ..StructuralModel: AbstractStructure, StructuralMaterials, StructuralBoundaryConditions


"""
An `Structure` object facilitates the process of assembling and creating the structural analysis. 
### Fields:
- `mesh`      -- Stores the structural mesh. 
- `materials` -- Stores the structural materials of the structure. 
- `elements`  -- Stores the structural elements of the structure.
- `bcs`       -- Stores the structural boundary conditions of the structure.
- `free_dofs` -- Stores the free degrees of freedom.
"""
mutable struct Structure{dim,MESH,MAT,E,NB,LB} <: AbstractStructure{dim,MAT,E}
    mesh::MESH
    materials::StructuralMaterials{MAT,E}
    bcs::StructuralBoundaryConditions{NB,LB}
    free_dofs::Vector{Dof}
    function Structure(
        mesh::MESH,
        materials::StructuralMaterials{MAT,E},
        bcs::StructuralBoundaryConditions{NB,LB},
        free_dofs::Vector{Dof}
    ) where {dim,MESH<:AbstractMesh{dim},MAT,E,NB,LB}
        return new{dim,MESH,MAT,E,NB,LB}(mesh, materials, bcs, free_dofs)
    end
end

"Constructor with  `StructuralMaterials` `materials`,  `StructuralBoundaryConditions` `bcs` 
and `AbstractMesh` `mesh` seting fixed dofs with `FixedDofBoundaryCondition` defined in `bcs`"
function Structure(
    mesh::AbstractMesh{dim},
    materials::StructuralMaterials{M,E},
    bcs::StructuralBoundaryConditions{NB,LB},
) where {dim,M,E,NB,LB}

    default_free_dofs = Vector{Dof}()
    for node_dofs in dofs(mesh)
        [push!(default_free_dofs, vec_dof...) for vec_dof in collect(values(node_dofs))]
    end

    fixed_dofs = compute_fixed_dofs(bcs, fixed_dof_bcs(bcs))

    deleteat!(default_free_dofs, findall(x -> x in fixed_dofs, default_free_dofs))

    return Structure(mesh, materials, bcs, default_free_dofs)
end

#### BoundaryConditions Applied to the structure
"Computes `Dof`s to delete given a `FixedDofBoundaryCondition` and a set of `StructuralBoundaryConditions` `bcs`."
function compute_fixed_dofs(bcs::StructuralBoundaryConditions, fbc::FixedDofBoundaryCondition)

    # Extract dofs to apply the bc
    fbc_dofs_symbols = dofs(fbc)

    # Extract nodes, faces and elements 
    entities = bcs[fbc]

    dofs_to_delete = Dof[]

    for dof_symbol in fbc_dofs_symbols
        dofs_entities = getindex.(dofs.(entities), dof_symbol)
        for component in components(fbc)
            dofs_to_delete_fbc = getindex.(dofs_entities, component)
            push!(dofs_to_delete, dofs_to_delete_fbc...)
        end
    end
    return dofs_to_delete
end

"Applies a `Vector` of `FixedDofBoundaryCondition` `f_bcs` and a set of `StructuralBoundaryConditions` `bcs`."
function compute_fixed_dofs(bcs::StructuralBoundaryConditions, f_bcs::Vector{<:FixedDofBoundaryCondition})
    dofs_to_delete = Dof[]
    [push!(dofs_to_delete, compute_fixed_dofs(bcs, fbc)...) for fbc in f_bcs]
    return dofs_to_delete
end
