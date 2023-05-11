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
struct PointsInterpolator{N<:AbstractNode,T<:Real,E<:AbstractElement,VE<:AbstractVector{E}} <:
       AbstractInterapolator
    "`Dictionary` with `Node`s as keys and the corresponding weights as values."
    node_to_weights::Vector{Dictionary{N,T}}
    "`Element` where the point is located."
    points_to_element::VE
end

"Return a `Vector` that maps the weight corresponding to each `Node`."
node_to_weights(points_interpolator::PointsInterpolator) = points_interpolator.node_to_weights

"Return a `Vector` that maps the `Element` where each `Point` is located."
points_to_element(points_interpolator::PointsInterpolator) = points_interpolator.points_to_element

"""
A `PointEvalHandler` facilitates the process of evaluating a solution at a given vector of points 
obtained at the `Node`s `Dof`s in a `Mesh`.
"""
struct PointEvalHandler{dim,T,PT<:Point{dim,T},VPT<:AbstractVector{PT},WT<:AbstractVector{T},
                        M<:AbstractMesh,
                        I<:AbstractInterapolator}
    "`Mesh` where the solution is obtained."
    mesh::M
    "Vector of test points."
    vec_points::Vector{PT}
    "Vector of indices of `vec_points` that lie inside the mesh."
    in_mesh_points_idx::Vector{Int64}
    "Vector of indices of `vec_points` that don't line inside the mesh."
    not_in_mesh_points_idx::Vector{Int64}
    "`Vector` of `Point`s to evaluate the solution."
    in_mesh_points::VPT
    "`Vector` of `Point`s outside the mesh."
    not_in_mesh_points::VPT
    "Vector of indices of `elements` of the mesh corresponding to each entry of `in_mesh_points_idx`."
    in_mesh_elements_idx::Vector{Int64}
    "Vector of weights per pair of `(element, point)` in the handler."
    weights::Vector{WT}
    "`Interpolator` object used to evaluate the solution."
    interpolator::I
end
function PointEvalHandler(mesh::AbstractMesh, point::P) where {T,P<:Point{T}}
    PointEvalHandler(mesh, [point])
end
function PointEvalHandler(mesh::AbstractMesh, vec_points::Vector{Vector{T}}) where {T}
    PointEvalHandler(mesh, [Point(p...) for p in vec_points])
end

"Return the vector of points where the solution will be evaluated."
points(peh::PointEvalHandler) = peh.in_mesh_points

"Return the vector of points that are not in the mesh."
not_in_mesh_points(peh::PointEvalHandler) = peh.not_in_mesh_points

"Return the `Mesh` where the solution is obtained."
mesh(peh::PointEvalHandler) = peh.mesh

"Return the `Interpolator` used to evaluate the solution."
interpolator(peh::PointEvalHandler) = peh.interpolator

abstract type SearchAlgorithm end

struct Serial <: SearchAlgorithm end
struct Threaded <: SearchAlgorithm end

"Constructor of a `PointEvalHandler` given a mesh and an array of points."
function PointEvalHandler(mesh::AbstractMesh, vec_points::Vector{PT};
                          alg::SearchAlgorithm=Threaded()) where {dim,T,PT<:Point{dim,T}}
    @assert dim ≤ 3 "Points must be 1D, 2D or 3D"

    # For each point, obtain the element(s) that it belongs to.
    # If a point belongs to more than one element, keep only the first matching element.
    # This is valid since the interpolation result will be the same for both elements.
    # This case occurs when the point is located on the boundary of two elements, face or node.
    in_mesh_points_idx, in_mesh_elements_idx = evaluate_points_in_mesh(mesh, vec_points, alg)

    # For each element found, compute the associated weight used for interpolation.
    # Assume that the mesh has a unique element type.
    nnodes = num_nodes(element(mesh, 1))
    weights = Vector{SVector{nnodes,T}}()
    @inbounds for (elem_idx, point_idx) in zip(in_mesh_elements_idx, in_mesh_points_idx)
        w = Elements.weights(element(mesh, elem_idx), vec_points[point_idx])
        push!(weights, w)
    end
    PointEvalHandler(mesh, vec_points, in_mesh_points_idx, in_mesh_elements_idx, weights)
end
function PointEvalHandler(mesh::AbstractMesh, vec_points::Vector{PT},
                          in_mesh_points_idx::Vector{Int64},
                          in_mesh_elements_idx::Vector{Int64},
                          weights::Vector{WT}) where {dim,T,PT<:Point{dim,T},WT<:AbstractVector{T}}
    # Subset of `vec_points` that belong to the mesh.
    in_mesh_points = view(vec_points, in_mesh_points_idx)

    # Subset of `vec_points` that do not belong to the mesh.
    not_in_mesh_points_idx = collect(1:length(vec_points))
    setdiff!(not_in_mesh_points_idx, in_mesh_points_idx)
    not_in_mesh_points = view(vec_points, not_in_mesh_points_idx)

    # Element where each point is located. At this stage, only one element per point is provided.
    points_to_element = view(elements(mesh), in_mesh_elements_idx)

    # Dictionary with nodes as keys and corresponding weights as values.
    num_in_mesh_points = length(in_mesh_elements_idx)
    node_to_weights = Vector{Dictionary{AbstractNode,T}}(undef, num_in_mesh_points)
    @inbounds for (i, elem_idx) in enumerate(in_mesh_elements_idx)
        elem = element(mesh, elem_idx)
        dict = dictionary([n => weights[i][j] for (j, n) in enumerate(nodes(elem))])
        node_to_weights[i] = dict
    end
    interpolator = PointsInterpolator(node_to_weights, points_to_element)

    PointEvalHandler(mesh, vec_points, in_mesh_points_idx, not_in_mesh_points_idx, in_mesh_points,
                     not_in_mesh_points, in_mesh_elements_idx, weights, interpolator)
end

function evaluate_points_in_mesh(mesh::AbstractMesh, vec_points::Vector{PT},
                                 ::Serial) where {dim,T,PT<:Point{dim,T}}
    in_mesh_points_idx = Vector{Int64}()
    in_mesh_elements_idx = Vector{Int64}()
    for (point_idx, point) in enumerate(vec_points)
        for (elem_idx, elem) in enumerate(elements(mesh))
            if point ∈ elem
                push!(in_mesh_points_idx, point_idx)
                push!(in_mesh_elements_idx, elem_idx)
                # Continue with next point, without checking the following elements.
                break
            end
        end
    end
    in_mesh_points_idx, in_mesh_elements_idx
end

function evaluate_points_in_mesh(mesh::AbstractMesh, vec_points::Vector{PT},
                                 ::Threaded) where {dim,T,PT<:Point{dim,T}}
    numthreads = Threads.nthreads()
    in_mesh_points_idx = [Vector{Int64}() for _ in 1:numthreads]
    in_mesh_elements_idx = [Vector{Int64}() for _ in 1:numthreads]
    Threads.@threads for point_idx in 1:length(vec_points)
        id = Threads.threadid()
        point = vec_points[point_idx]
        for (elem_idx, elem) in enumerate(elements(mesh))
            if point ∈ elem
                push!(in_mesh_points_idx[id], point_idx)
                push!(in_mesh_elements_idx[id], elem_idx)
                # Continue with next point, without checking the following elements.
                break
            end
        end
    end
    return reduce(vcat, in_mesh_points_idx), reduce(vcat, in_mesh_elements_idx)
end

end # module
