"""
Module to export ONSAS data structures to `.vtk` files suitable for [Paraview](https://www.paraview.org/).
"""
module VTK

using Reexport
using WriteVTK

using ..Meshes
using ..Nodes
using ..Utils
using ..Entities
using ..Tetrahedrons
using ..Solutions
using ..Structures
using ..StructuralAnalyses

export VTKMeshFile, create_vtk_grid, write_vtk

struct VTKMeshFile{VTK<:WriteVTK.DatasetFile}
    "`WriteVTK` `VTK` native file."
    vtk::VTK
end
function VTKMeshFile(filename::String, Mesh::AbstractMesh; kwargs...)
    vtk = create_vtk_grid(filename, Mesh; kwargs...)
    VTKMeshFile(vtk)
end
# Functor handler allowing do syntax
function VTKMeshFile(f::Function, args...; kwargs...)
    vtk = VTKMeshFile(args...; kwargs...)
    try
        f(vtk)
    finally
        close(vtk)
    end
    vtk
end
function Base.show(io::IO, ::MIME"text/plain", (; vtk)::VTKMeshFile)
    open_str = isopen(vtk) ? "open" : "closed"
    filename = vtk.path
    ncells = vtk.Ncls
    nnodes = vtk.Npts
    print(io,
          "• VTKMeshFile file \"$(filename)\" is $open_str with $nnodes nodes and $ncells cells.")
end

function create_vtk_grid(filename::String, mesh::AbstractMesh{dim}; kwargs...) where {dim}
    cls = WriteVTK.MeshCell[]
    for elem in elements(mesh)
        VTK_elem_type = to_vtkcell_type(elem)
        push!(cls, WriteVTK.MeshCell(VTK_elem_type, to_vtk_cell_nodes(elem, mesh)))
    end
    coords = node_matrix(mesh)
    WriteVTK.vtk_grid(filename, coords, cls; kwargs...)
end

function Base.close(vtk::VTKMeshFile)
    WriteVTK.vtk_save(vtk.vtk)
end

to_vtkcell_type(::Tetrahedron) = VTKCellTypes.VTK_TETRA
function to_vtk_cell_nodes(e::AbstractElement, msh::AbstractMesh)
    connectivity(msh)[findfirst(==(e), elements(msh))]
end

function to_vtk!(v, x::Vector{dim}) where {dim}
    v[1:dim] .= x
end
function to_vtk!(v, x::AbstractMatrix)
    v[1:length(vog)] .= vogit(x)
    v
end

function _vtk_write_node_data(vtk::WriteVTK.DatasetFile,
                              nodedata::Vector{<:Real},
                              name::AbstractString)
    WriteVTK.vtk_point_data(vtk, nodedata, name)
end
function _vtk_write_node_data(vtk::WriteVTK.DatasetFile,
                              nodedata::Matrix{<:Real},
                              name::AbstractString;
                              component_names=nothing)
    WriteVTK.vtk_point_data(vtk, nodedata, name; component_names=component_names)
end

"""
Generate a VTK file given a solution struct.
Currently only `displacements` are exported for each time point.
"""
function write_vtk(sol::AbstractSolution, filename::String)
    msh = mesh(structure(analysis(sol)))
    nodes_mat = node_matrix(msh)
    connec_mat = connectivity(msh)
    num_elem = num_elements(msh)
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, connec_mat[e]) for e in 1:num_elem]

    n_times = length(displacements(sol))
    n_dofs = num_dofs(msh)
    mypad = Integer(ceil(log10(n_times))) + 1
    for i in 1:n_times
        pdata = displacements(sol)[i]
        nodes = nodes_mat + reshape(pdata, (3, Integer(n_dofs / 3))) # TO DO generalize for angle dof cases
        filename_i = filename * string(i; pad=mypad)
        WriteVTK.vtk_grid(filename_i, nodes, cells) do vtk
            vtk["Displacements", VTKPointData()] = pdata
        end
    end
    @info "VTK output written to $filename"
    filename
end

end
