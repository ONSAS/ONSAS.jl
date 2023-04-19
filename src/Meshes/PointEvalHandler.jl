using Reexport: @reexport

using Dictionaries: Dictionary
using ..Meshes: AbstractMesh, elements
using ..Elements: AbstractNode, AbstractElement, coordinates, index, weights

export PointEvalHandler, points, interpolator, mesh
export PointsInterpolator, node_to_weights, points_to_element

""" PointEvalHandler struct.
A `PointEvalHandler` facilitates the process of evaluating a solution at a given vector of points 
obtained at the `Node`s `Dof`s in a `Mesh`.
### Fields:
- `vec_points`   -- stores the points where the solution will be evaluated.
- `mesh`         -- stores the mesh where the solution is obtained.
- `interpolator` -- stores the interpolator used to evaluate the solution.
"""
struct PointEvalHandler{I,M,P,VP<:Vector{P}}
    vec_points::VP
    mesh::M
    interpolator::I
end

"Returns the `Vector` of `Point`s where the solution will be evaluated."
points(peh::PointEvalHandler) = peh.vec_points

"Returns the `Mesh` where the solution is obtained."
mesh(peh::PointEvalHandler) = peh.mesh

"Returns the `Interpolator` used to evaluate the solution."
interpolator(peh::PointEvalHandler) = peh.interpolator

"""PointsInterpolator struct.
A `PointsInterpolator` struct stores the weights nodes and elements needed to interpolate the solution at a given point.
The index of each `Vector` is the index in the `Vector` of `Point`s.
### Fields:
- `node_to_weights` -- stores a `Dictionary` with `Node`s as keys and the corresponding weights as values.
- `point_to_element` -- stores the `Element` where the point is located.
"""
struct PointsInterpolator{N<:AbstractNode,T,E<:AbstractElement}
    node_to_weights::Vector{Dictionary{N,T}}
    points_to_element::Vector{E}
end

"Returns a `Vector` that maps the weight corresponding to each `Node`."
node_to_weights(points_interpolator::PointsInterpolator) = points_interpolator.node_to_weights

"Returns a `Vector` that maps the `Element` where each `Point` is located."
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
            node_to_weights_p = dictionary([node => weights_p[node_index] for (node_index, node) in enumerate(e_nodes)])
            point_interpolators[point_index] = node_to_weights_p

            # Add the point to the interpolated points
            push!(interpolated_vec_points_indexes, point_index)

            # Add the element
            point_to_element[point_index] = e
        end

    end

    @assert length(interpolated_vec_points_indexes) == length(vec_points) "
    There are $(length(vec_points) - length(interpolated_vec_points_indexes)) points outside the mesh .
    "

    return PointEvalHandler(vec_points, mesh, PointsInterpolator(point_interpolators, point_to_element))
end

"Constructor of a `PointEvalHandler` from a `Mesh` and a `Point` ."
PointEvalHandler(mesh::AbstractMesh, point::P) where {T,P<:Point{T}} = PointEvalHandler(mesh, [point])
