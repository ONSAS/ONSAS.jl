using Reexport: @reexport

using Dictionaries: Dictionary
using ..Meshes: AbstractMesh, elements
using ..Elements: coordinates, index, weights

const Point{T} = Union{<:AbstractVector{T},<:NTuple{dim,T}} where {dim,T<:Real}

export PointEvalHandler, points, interpolator, mesh


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

"Constructor of a `PointEvalHandler` from a `Mesh` , the `Dof` symbols into a `Vector`and a `AbstractVector` of `Point`s . "
function PointEvalHandler(mesh::AbstractMesh, vec_points::AbstractVector{P}) where {T,P<:Point{T}}

    @assert length(first(vec_points)) ≤ 3 "Points must be 1D, 2D or 3D"

    # Init a list of Dictionaries to store the weights and nodes needed to interpolate the file at each point
    point_interpolators = Vector{Dictionary{AbstractNode,Float64}}(undef, length(vec_points))

    # Init a list to store the indexes of the points that have been already interpolated
    interpolated_vec_points_indexes = Int[]
    sizehint!(interpolated_vec_points_indexes, length(vec_points))

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
        end

    end

    @assert length(interpolated_vec_points_indexes) == length(vec_points) "There are points outside the mesh."

    return PointEvalHandler(vec_points, mesh, point_interpolators)
end

