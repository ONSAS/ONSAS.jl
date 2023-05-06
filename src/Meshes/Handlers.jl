module Handlers

using Reexport, Dictionaries

using ..Meshes
using ..Elements

export PointEvalHandler, points, not_in_mesh_points, interpolator, mesh, PointsInterpolator,
       node_to_weights, points_to_element

"""
A `PointEvalHandler` facilitates the process of evaluating a solution at a given vector of points 
obtained at the `Node`s `Dof`s in a `Mesh`.
"""
struct PointEvalHandler{I,M<:AbstractMesh,P,VP1<:AbstractVector{P},VP2<:AbstractVector{P}}
    "`Vector` of `Point`s to evaluate the solution."
    in_mesh_points::VP1
    "`Vector` of `Point`s outside the mesh."
    not_in_mesh_points::VP2
    "`Mesh` where the solution is obtained."
    mesh::M
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

"""
A `PointsInterpolator` struct stores the weights nodes and elements needed 
    to interpolate the solution at a given point.The index of each `Vector` 
    is the index in the `Vector` of `Point`s.
"""
struct PointsInterpolator{N<:AbstractNode,T,E<:AbstractElement}
    "`Dictionary` with `Node`s as keys and the corresponding weights as values."
    node_to_weights::Vector{Dictionary{N,T}}
    "`Element` where the point is located."
    points_to_element::Vector{E}
end

"Return a `Vector` that maps the weight corresponding to each `Node`."
node_to_weights(points_interpolator::PointsInterpolator) = points_interpolator.node_to_weights

"Return a `Vector` that maps the `Element` where each `Point` is located."
points_to_element(points_interpolator::PointsInterpolator) = points_interpolator.points_to_element

"Constructor of a `PointEvalHandler` from a `Mesh` and a `AbstractVector` of `Point`s ."
function PointEvalHandler(mesh::AbstractMesh, vec_points::AbstractVector{P}) where {T,P<:Point{T}}
    @assert length(first(vec_points)) ≤ 3 "Points must be 1D, 2D or 3D"

    # Init a list of Dictionaries to store the weights and nodes needed to interpolate the file at each point
    point_interpolators = Vector{Dictionary{AbstractNode,Float64}}(undef, length(vec_points))

    # Init a list to store the indexes of the points that have been already interpolated
    interpolated_vec_points_indexes = Int[]
    sizehint!(interpolated_vec_points_indexes, length(vec_points))

    # Init a list to store the elements where the points are located
    point_to_element = Vector{AbstractElement}(undef, length(vec_points))

    # Check for every element which points are inside the element and compute the 
    # weights using the shape function of the element 
    for e in elements(mesh)
        e_nodes = nodes(e)

        # If the points have been already evaluated break out
        length(interpolated_vec_points_indexes) == length(vec_points) && break

        for (point_index, p) in enumerate(vec_points)

            # If the point has been already evaluated skip it
            point_index ∈ interpolated_vec_points_indexes && continue
            # If the point is not inside the element skip it
            p ∉ e && continue

            weights_p = weights(e, p)

            # Store the weights for each node
            node_to_weights_p = dictionary([node => weights_p[node_index]
                                            for (node_index, node) in enumerate(e_nodes)])
            point_interpolators[point_index] = node_to_weights_p

            # Add the point to the interpolated points
            push!(interpolated_vec_points_indexes, point_index)

            # Add the element
            point_to_element[point_index] = e
        end
    end

    if length(interpolated_vec_points_indexes) == length(vec_points)
        PointEvalHandler(vec_points, Vector{P}(undef, 0), mesh,
                         PointsInterpolator(point_interpolators, point_to_element))
    else
        @warn "There are $(length(vec_points) - length(interpolated_vec_points_indexes)) points outside the mesh ."
        # Pick in mesh and no in mesh points
        sort!(interpolated_vec_points_indexes)
        in_mesh_points = view(vec_points, interpolated_vec_points_indexes)
        not_in_mesh_indexes = deleteat!(collect(1:length(vec_points)),
                                        interpolated_vec_points_indexes)
        not_in_mesh_points = view(vec_points, not_in_mesh_indexes)

        # Remove undef in interpolator object
        deleteat!(point_interpolators, not_in_mesh_indexes)
        deleteat!(point_to_element, not_in_mesh_indexes)
        interpolator = PointsInterpolator(point_interpolators, point_to_element)
        # Return the PointEvalHandler
        PointEvalHandler(in_mesh_points, not_in_mesh_points, mesh, interpolator)
    end
end

"Constructor of a `PointEvalHandler` from a `Mesh` and a `Point` ."
function PointEvalHandler(mesh::AbstractMesh, point::P) where {T,P<:Point{T}}
    PointEvalHandler(mesh, [point])
end
function PointEvalHandler(mesh::AbstractMesh, vec_points::Vector{Vector{T}}) where {T}
    PointEvalHandler(mesh, [Point(p...) for p in vec_points])
end

end # module
