using ..Meshes: AbstractMesh, elements
using ..Elements: coordinates, index, weights

using Dictionaries: Dictionary
import SparseArrays: sparse

const Point{T} = Union{<:AbstractVector{T},<:NTuple{dim,T}} where {dim,T<:Real}

export PointEvalHandler, node_to_elements


""" PointEvalHandler struct.
A `PointEvalHandler` facilitates the process of evaluating a solution at a given vector of points 
obtained at the `Node`s `Dof`s in a `Mesh`.
### Fields:
- `vec_points`   -- stores the points where the solution will be evaluated.
- `mesh`         -- stores the mesh where the solution is obtained.
- `interpolator` -- stores the interpolator used to evaluate the solution.
- `dof_symbols`  -- stores the `Dof` symbols that will be used to interpolate the solution.
"""
struct PointEvalHandler{I,M,P,VP<:Vector{P}}
    vec_points::VP
    mesh::M
    interpolator::I
    dof_symbols::Vector{Symbol}
end

"Returns the `Vector` of `Point`s where the solution will be evaluated."
vec_points(peh::PointEvalHandler) = peh.vec_points

"Returns the `Mesh` where the solution is obtained."
mesh(peh::PointEvalHandler) = peh.mesh

"Returns the `Interpolator` used to evaluate the solution."
interpolator(peh::PointEvalHandler) = peh.interpolator

"Returns the interpolator matrix for the dof symbol `dof_symbol`."
Base.getindex(peh::PointEvalHandler, dof_symbol::Symbol) = peh.interpolator[dof_symbol]

"Returns the `Dof` symbols that will be used to interpolate the solution."
dof_symbols(peh::PointEvalHandler) = peh.dof_symbols




""" SparseInterpolator struct.
A `SparseInterpolator` helps to build an sparse interpolator matrix.
### Fields:
- `dofs_indexes`    -- stores dof indexes that will be used to interpolate the solution.
- `points_indexes`  -- stores the points indexes where the solution will be evaluated.
- `weights`         -- stores the interpolation weights.
"""
struct SparseInterpolatorAssembler{T}
    dofs_indexes::Vector{Int}
    points_indexes::Vector{Int}
    weights::Vector{T}
end

"Constructor for `CSCSparseMatrix` from a `SparseInterpolator` ."
sparse(sia::SparseInterpolatorAssembler, num_dofs_mesh::Int, num_points::Int) =
    sparse(sia.points_indexes, sia.dofs_indexes, sia.weights, num_points, num_dofs_mesh, +)

"Constructor of an `SparseInterpolatorAssembler` with size of the sparse matrix `N` ."
function SparseInterpolatorAssembler(N::Integer)
    dofs_indexes = Int[]
    nodes_to_interpolate = Int[]
    weights = Float64[]

    sizehint!(dofs_indexes, N)
    sizehint!(nodes_to_interpolate, N)
    sizehint!(weights, N)

    SparseInterpolatorAssembler(dofs_indexes, nodes_to_interpolate, weights)
end

"Constructor of a `PointEvalHandler` from a `Mesh` , the `Dof` symbols into a `Vector`and a `AbstractVector` of `Point`s . "
function PointEvalHandler(
    mesh::AbstractMesh, dof_symbols::Vector{Symbol}, vec_points::AbstractVector{P}
) where {T,P<:Point{T}}

    @assert length(first(vec_points)) ≤ 3 "Points must be 1D, 2D or 3D"

    num_dofs_mesh = num_dofs(mesh)
    num_points = length(vec_points)

    interpolated_vec_points_indexes = Vector{Int}()
    sizehint!(interpolated_vec_points_indexes, num_points)

    # Stores an interpolation matrix for each degree of symbol
    interpolation_assembler = dictionary(
        [dof_symbol => SparseInterpolatorAssembler(num_dofs(mesh)) for dof_symbol in dof_symbols]
    )

    # Check for every element which points are inside the element, interpolate 
    # and store indexes and weights
    for e in elements(mesh)

        # If the points have been already evaluated break out
        length(interpolated_vec_points_indexes) == length(vec_points) && break

        element_coords = coordinates.(nodes(e))
        min = first(minimum(element_coords, dims=1))
        max = first(maximum(element_coords, dims=1))

        for (point_index, p) in enumerate(vec_points)
            # If point is outside the box check and the element check the next point
            all(p .>= min) && all(p .<= max) || continue

            # Check if the point is inside the element
            # p ⊄ e || continue

            # If inside compute weights
            weights_p = weights(e, p)

            # Build a dictionary that sotres as keys the dof symbols and as values the assembler object to build a sparse 
            # matrix (dof_index, point_index) -> weight
            for dof_symbol in dof_symbols
                # Get the assembler object for the dof_symbol
                dof_assembler = interpolation_assembler[dof_symbol]

                for (node_element_index, node) in enumerate(nodes(e))

                    node_weight = weights_p[node_element_index]
                    node_dofs = dofs(node)[dof_symbol]
                    num_dofs_node = length(node_dofs)
                    # Add the weights and dofs that are required to interpolate
                    # the solution in p
                    push!(dof_assembler.dofs_indexes, index.(node_dofs)...)
                    # Add the point index (in vec_points) that is being interpolated 
                    push!(dof_assembler.points_indexes, repeat([point_index], num_dofs_node)...)
                    # Add the weight to the node
                    push!(dof_assembler.weights, repeat([node_weight], num_dofs_node)...)

                end

            end

            # Add the point to the interpolated points
            push!(interpolated_vec_points_indexes, point_index)
        end

    end
    interpolator = dictionary(
        [dof_symbol => sparse(interpolation_assembler[dof_symbol], num_dofs_mesh, num_points) for dof_symbol in dof_symbols]
    )
    return PointEvalHandler(vec_points, mesh, interpolator, dof_symbols)
end

