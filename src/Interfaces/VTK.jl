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

export VTKMeshFile, create_vtk_grid, write_node_data, write_vtk

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

function _vtk_write_node_data(vtk::WriteVTK.DatasetFile,
                              nodal_data::Vector{<:Real},
                              name::AbstractString)
    WriteVTK.vtk_point_data(vtk, nodal_data, name)
end
function _vtk_write_node_data(vtk::WriteVTK.DatasetFile,
                              nodal_data::Matrix{<:Real},
                              name::AbstractString;
                              component_names=nothing)
    WriteVTK.vtk_point_data(vtk, nodal_data, name; component_names=component_names)
end
function write_node_data(vtk::VTKMeshFile, nodedata, name; kwargs...)
    _vtk_write_node_data(vtk.vtk, nodedata, name; kwargs...)
    vtk
end

to_vtk(x::Vector) = x
VTX_VOIGT_ORDER = [1, 6, 5, 9, 2, 4, 8, 7, 3]
to_vtk(x::AbstractMatrix) = view(x, VTX_VOIGT_ORDER)

function default_dof_fields(sol::AbstractSolution)
    ns = nodes(mesh(structure(analysis(sol))))
    nf = collect(keys(dofs(first(ns))))
    push!(nf, :σ)
    push!(nf, :ϵ)
    nf
end

const POINT_FIELDS = [:u, :θ]
const CELL_FIELDS = Dict(:σ => (3, 3), :ϵ => (3, 3))
const FIELD_NAMES = Dict(:u => "Displacement", :θ => "Rotation", :σ => "Stress", :ϵ => "Strain")

"""
Generate a VTK file given a solution struct.
Currently only `displacements` are exported for each time point.
"""
function write_vtk(sol::AbstractSolution,
                   filename::String;
                   vf::Vector{Field}=default_dof_fields(sol),
                   time_index::Integer)
    msh = mesh(structure(analysis(sol)))

    nodal_data = Dict(field => zeros(num_dofs(msh, field))
                      for field in vf if field ∈ POINT_FIELDS)
    for node in nodes(msh)
        for field in vf
            if field ∈ POINT_FIELDS
                ndofs = dofs(node, field)
                node_dof_data = displacements(sol, ndofs, time_index)
                nodal_data[field][ndofs] .= to_vtk(node_dof_data)
            end
        end
    end

    cell_data = Dict(field => zeros(prod(CELL_FIELDS[field]) * num_elements(msh))
                     for field in vf if field ∈ keys(CELL_FIELDS))

    for (i, elem) in enumerate(elements(msh))
        for field in vf
            if haskey(CELL_FIELDS, field)
                num_field_comps = prod(CELL_FIELDS[field])
                if field ∈ keys(CELL_FIELDS)
                    start = (i - 1) * num_field_comps + 1
                    finish = start + num_field_comps - 1
                    if field == :σ
                        σ = stress(sol, elem, time_index)
                        cell_data[field][start:finish] .= to_vtk(σ)
                    elseif field == :ϵ
                        ϵ = strain(sol, elem, time_index)
                        cell_data[field][start:finish] .= to_vtk(ϵ)
                    end
                end
            end
        end
    end

    VTKMeshFile(filename, msh) do vtx
        for (field, data) in nodal_data
            field_name = FIELD_NAMES[field]
            write_node_data(vtx, data, field_name)
        end
        for (field, data) in cell_data
            field_name = FIELD_NAMES[field]
            write_cell_data(vtx, data, field_name)
        end
    end
    @info "VTK output written to $filename"
    filename
end

end
