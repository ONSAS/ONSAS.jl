module Handlers

using Reexport, Dictionaries
using StaticArrays: SVector

using ..Meshes
using ..Elements

export PointEvalHandler, points, not_in_mesh_points, interpolator, mesh, PointsInterpolator,
       node_to_weights, points_to_element

abstract type AbstractInterapolator end

"""
A `PointsInterpolator` struct stores the weights nodes and elements needed
to interpolate the solution at a given point.The index of each `Vector`
is the index in the `Vector` of `Point`s.
"""
struct PointsInterpolator{N<:AbstractNode,T,E<:AbstractElement} <: AbstractInterapolator
    "`Dictionary` with `Node`s as keys and the corresponding weights as values."
    node_to_weights::Vector{Dictionary{N,T}}
    "`Element` where the point is located."
    points_to_element::Vector{E}
end

"Return a `Vector` that maps the weight corresponding to each `Node`."
node_to_weights(points_interpolator::PointsInterpolator) = points_interpolator.node_to_weights

"Return a `Vector` that maps the `Element` where each `Point` is located."
points_to_element(points_interpolator::PointsInterpolator) = points_interpolator.points_to_element

"""
A `PointEvalHandler` facilitates the process of evaluating a solution at a given vector of points 
obtained at the `Node`s `Dof`s in a `Mesh`.
"""
struct PointEvalHandler{dim,T,PT<:Point{dim,T},WT<:AbstractVector{T},M<:AbstractMesh,
                        I<:AbstractInterapolator}
    "`Mesh` where the solution is obtained."
    mesh::M
    "Vector of test points."
    vec_points::Vector{PT}
    "Vector of indices of `vec_points` that lie inside the mesh."
    in_mesh_points_idx::Vector{Int64}
    "Vector of indices of `elements` of the mesh corresponding to each entry of `in_mesh_points_idx`."
    in_mesh_elements_idx::Vector{Int64}
    "Vector of weights per pair of `(element, point)` in the handler."
    weights::Vector{WT}
    "`Vector` of `Point`s to evaluate the solution."
    in_mesh_points::Vector{PT}
    "`Vector` of `Point`s outside the mesh."
    not_in_mesh_points::Vector{PT}
    "`Interpolator` object used to evaluate the solution."
    interpolator::I
end

"Return the `Vector` of `Point`s where the solution will be evaluated."
points(peh::PointEvalHandler) = peh.in_mesh_points

"Return the `Vector` of `Point`s that are not in the mesh."
not_in_mesh_points(peh::PointEvalHandler) = peh.not_in_mesh_points

"Return the `Mesh` where the solution is obtained."
mesh(peh::PointEvalHandler) = peh.mesh

"Return the `Interpolator` used to evaluate the solution."
interpolator(peh::PointEvalHandler) = peh.interpolator

"Constructor of a `PointEvalHandler` given a mesh and an array of points."
function PointEvalHandler(mesh::AbstractMesh, vec_points::Vector{PT}) where {dim,T,PT<:Point{dim,T}}
    @assert dim ≤ 3 "Points must be 1D, 2D or 3D"

    # For each point, obtain the elements that it belongs to.
    in_mesh_points_idx = Vector{Int64}()
    in_mesh_elements_idx = Vector{Int64}()
    for (point_idx, point) in enumerate(vec_points)
        for (elem_idx, elem) in enumerate(elements(mesh))
            if point ∈ elem
                push!(in_mesh_points_idx, point_idx)
                push!(in_mesh_elements_idx, elem_idx)
                # TODO Add `in_mesh_points` in this loop if `point_idx` doesn't belong to `in_mesh_points_idx`. 
            end
        end
    end

    # For each element found, compute the associated weight used for interpolation.
    n = length(in_mesh_points_idx)
    # Assuming tetrahedron mesh (4 nodes per element); may otherwise use inner `Vector` if needed.
    weights = Vector{SVector{4,T}}()
    @inbounds for (elem_idx, point_idx) in zip(in_mesh_elements_idx, in_mesh_points_idx)
        w = Elements.weights(element(mesh, elem_idx), vec_points[point_idx])
        push!(weights, w)
    end
    PointEvalHandler(mesh, vec_points, in_mesh_points_idx, in_mesh_elements_idx, weights)
end
function PointEvalHandler(mesh::AbstractMesh, point::P) where {T,P<:Point{T}}
    PointEvalHandler(mesh, [point])
end
function PointEvalHandler(mesh::AbstractMesh, vec_points::Vector{Vector{T}}) where {T}
    PointEvalHandler(mesh, [Point(p...) for p in vec_points])
end

# TODO Cleanup redundant fields.
function PointEvalHandler(mesh::AbstractMesh, vec_points::Vector{PT},
                          in_mesh_points_idx::Vector{Int64},
                          in_mesh_elements_idx::Vector{Int64},
                          weights::Vector{WT}) where {dim,T,PT<:Point{dim,T},WT<:AbstractVector{T}}
    # Subset of `vec_points` that belong to the mesh.
    unique_in_mesh_points_idx = unique(in_mesh_points_idx)
    num_in_points = length(unique_in_mesh_points_idx)
    in_mesh_points = Vector{PT}(undef, num_in_points)
    for (i, idx) in enumerate(unique_in_mesh_points_idx)
        in_mesh_points[i] = vec_points[idx]
    end
    # unique_in_mesh_points_idx = unique(in_mesh_points_idx)
    # in_mesh_points = vec_points[unique_in_mesh_points_idx]
    # setdiff!(...)

    # Subset of `vec_points` that do not belong to the mesh.
    num_not_in_mesh_points = length(vec_points) - num_in_points
    not_in_mesh_points = Vector{PT}()
    sizehint!(not_in_mesh_points, num_not_in_mesh_points)
    for (i, p) in enumerate(vec_points)
        if i ∉ in_mesh_points_idx
            push!(not_in_mesh_points, p)
        end
    end

    # Main.@infiltrate
    # Element where each point is located.
    points_to_element = Vector{AbstractElement}()
    for (i, idx) in enumerate(in_mesh_elements_idx)
        push!(points_to_element, element(mesh, idx))
    end

    # Dictionary with nodes as keys and corresponding weights as values.
    node_to_weights = Vector{Dictionary{AbstractNode,T}}()
    sizehint!(node_to_weights, 4 * num_in_points)
    for (i, elem_idx) in enumerate(in_mesh_elements_idx)
        elem = element(mesh, elem_idx)
        @show nodes(elem)
        @show weights[i]
        dict = dictionary([n => weights[i][j] for (j, n) in enumerate(nodes(elem))])
        push!(node_to_weights, dict)
    end
    interpolator = PointsInterpolator(node_to_weights, points_to_element)

    PointEvalHandler(mesh, vec_points, in_mesh_points_idx, in_mesh_elements_idx, weights,
                     in_mesh_points, not_in_mesh_points, interpolator)
end

end # module
