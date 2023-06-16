"""
Module defining methods and types to handle mesh operations.
"""
module Handlers

using Reexport, Dictionaries, StaticArrays, LazySets

using ..Meshes
using ..Entities
using ..Nodes

@reexport import ..Utils: mesh

export AbstractHandler, PointEvalHandler, points, not_in_mesh_points, mesh, interpolator

""" 
An `AbstractHandler` object facilitates the process to handle meshes.
"""
abstract type AbstractHandler{dim,M<:AbstractMesh} end

"""
A `PointEvalHandler` facilitates the process of evaluating a solution 
at a given vector of points obtained at the `Node`'s `Dof` in the `AbstractMesh`.
"""
struct PointEvalHandler{dim,
                        T<:Real,
                        P<:Point{dim,T},
                        VP<:AbstractVector{P},
                        WT<:AbstractVector{T},
                        M<:AbstractMesh,
                        I<:AbstractInterpolator}
    "`Mesh` where the solution is obtained."
    mesh::M
    "Vector of test points."
    points::Vector{P}
    "Vector of indices of `points` that lie inside the mesh."
    in_mesh_points_idx::Vector{Int}
    "Vector of indices of `points` that don't line inside the mesh."
    not_in_mesh_points_idx::Vector{Int}
    "`Vector` of `Point`s to evaluate the solution."
    in_mesh_points::VP
    "`Vector` of `Point`s outside the mesh."
    not_in_mesh_points::VP
    "Vector of indices of `elements` of the mesh corresponding to each entry of `in_mesh_points_idx`."
    in_mesh_elements_idx::Vector{Int}
    "`Interpolator` object used to evaluate the solution."
    interpolator::I
    # Check the lengths
    function PointEvalHandler(mesh::M,
                              points::VP,
                              in_mesh_elements_idx::Vector{Int},
                              not_in_mesh_points_idx::Vector{Int},
                              in_mesh_points::VIP,
                              not_in_mesh_points::VNP,
                              in_mesh_elements_idx::Vector{Int},
                              interpolator::I) where {dim,
                                                      T<:Real,
                                                      P<:Point{dim,T},
                                                      VP<:AbstractVector{P},
                                                      M<:AbstractMesh,
                                                      I<:AbstractInterpolator}
        @assert length(points) == length(in_mesh_points) + length(not_in_mesh_points)
        new{dim,T,P,VP,WT,M,I}(mesh, points, in_mesh_points, not_in_mesh_points,
                               in_mesh_elements_idx, interpolator)
    end
end

"Constructor of a `PointEvalHandler` given a mesh and `Point`."
function PointEvalHandler(mesh::AbstractMesh, point::P) where {T,P<:Point{T}}
    PointEvalHandler(mesh, [point])
end

"Constructor of a `PointEvalHandler` given a mesh and `Vector` of `Point`s."
function PointEvalHandler(mesh::AbstractMesh, points::Vector{Vector{T}}) where {T}
    PointEvalHandler(mesh, [Point(p...) for p in points])
end

"Return the vector of points where the solution will be evaluated."
points(peh::PointEvalHandler) = peh.in_mesh_points

"Return the vector of points that are not in the mesh."
not_in_mesh_points(peh::PointEvalHandler) = peh.not_in_mesh_points

"Return the `Mesh` where the solution is obtained."
mesh(peh::PointEvalHandler) = peh.mesh

"Return the `Interpolator` used to evaluate the solution."
interpolator(peh::PointEvalHandler) = peh.interpolator

"Constructor of a `PointEvalHandler` given a `AbstractMesh` `mesh` and `Vector` of `Points`.
Moreover, the `SearchAlgorithm` `alg` can be specified to search the elements that contain the points."
function PointEvalHandler(mesh::AbstractMesh, points::Vector{P};
                          alg::AbstractSearchAlgorithm=Serial()) where {dim,T,P<:Point{dim,T}}
    @assert dim ≤ 3 "Points must be 1D, 2D or 3D"

    # For each point, obtain the element(s) that it belongs to.
    # If a point belongs to more than one element, keep only the first matching element.
    # This is valid since the interpolation result will be the same for both elements.
    # This case occurs when the point is located on the boundary of two elements, face or node.
    in_mesh_points_idx, in_mesh_elements_idx = evaluate_points_in_mesh(mesh, points, alg)

    # For each element found, compute the associated weight used for interpolation.
    # Assume that the mesh has a unique element type.
    nnodes = num_nodes(Meshes.element(mesh, 1))
    weights = Vector{SVector{nnodes,T}}(undef, length(in_mesh_elements_idx))
    @inbounds for (i, (elem_idx, point_idx)) in
                  enumerate(zip(in_mesh_elements_idx, in_mesh_points_idx))
        w = Elements.weights(Meshes.element(mesh, elem_idx), points[point_idx])
        weights[i] = w
   end
    PointEvalHandler(mesh, points, in_mesh_points_idx, in_mesh_elements_idx, weights)
end

"Constructor of a `PointEvalHandler` given a `AbstractMesh` `mesh` and `Vector` of `Points`.
Moreover, the `SearchAlgorithm` `alg` can be specified to search the elements that contain the points."
function PointEvalHandler(mesh::AbstractMesh, points::Vector{P},
                          in_mesh_points_idx::Vector{Int},
                          in_mesh_elements_idx::Vector{Int},
                          weights::Vector{WT}) where {dim,T,P<:Point{dim,T},WT<:AbstractVector{T}}
    # Subset of `points` that belong to the mesh.
    in_mesh_points = view(points, in_mesh_points_idx)

    # Subset of `points` that do not belong to the mesh.
    not_in_mesh_points_idx = collect(1:length(points))
    setdiff!(not_in_mesh_points_idx, in_mesh_points_idx)
    not_in_mesh_points = view(points, not_in_mesh_points_idx)

    # Element where each point is located. At this stage, only one element per point is provided.
    points_to_element = view(elements(mesh), in_mesh_elements_idx)

    # Dictionary with nodes as keys and corresponding weights as values.
    num_in_mesh_points = length(in_mesh_elements_idx)
    node_to_weights = Vector{Dictionary{AbstractNode{dim},T}}(undef, num_in_mesh_points)
    @inbounds for (i, elem_idx) in enumerate(in_mesh_elements_idx)
        elem = Meshes.element(mesh, elem_idx)
        dict = dictionary([n => weights[i][j] for (j, n) in enumerate(nodes(elem))])
        node_to_weights[i] = dict
    end
    interpolator = FEMInterpolator(view(points, in_mesh_points_idx),
                                   node_to_weights,
                                   points_to_element)
    Main.@infiltrate
    PointEvalHandler(mesh, points, in_mesh_points_idx, not_in_mesh_points_idx, in_mesh_points,
                     not_in_mesh_points, in_mesh_elements_idx, interpolator)
end

end # module
