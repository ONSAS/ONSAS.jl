using ..Meshes: AbstractMesh, elements
using ..Elements: coordinates, index

using Dictionaries: Dictionary
using SparseArrays: SparseMatrixCSC

const Point{T} = Vector{T} where {T<:Real}

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

""" SparseInterpolator struct.
A `SparseInterpolator` helps to build an sparse interpolator matrix.
### Fields:
- `dofs_indexes`          -- stores dof indexes that will be used to interpolate the solution.
- `nodes_to_interpolate`  -- stores the nodes indexes where the solution will be evaluated.
- `weights`               -- stores the interpolation weights.
"""
struct SparseInterpolatorAssembler{T}
    dofs_indexes::Vector{Int}
    nodes_to_interpolate::Vector{Int}
    weights::Vector{T}
end

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

function PointEvalHandler(
    mesh::AbstractMesh, dof_symbols::Vector{Symbol}, vec_points::Vector{P}
) where {T,P<:Point{T}}

    interpolated_vec_points = Vector{P}()
    sizehint!(interpolated_vec_points, length(vec_points))

    # Stores an interpolation matrix for each degree of symbol
    interpolation_assembler = dictionary(
        [dof_symbol => SparseInterpolatorAssembler(length(vec_points)) for dof_symbol in dof_symbols]
    )

    # Check for every element which points are inside the element, interpolate 
    # and store indexes and weights
    for e in elements(mesh)

        # If the points have been already evaluated break out
        length(interpolated_vec_points) == length(vec_points) && break

        element_coords = coordinates.(nodes(e))
        min = first(minimum(element_coords, dims=1))
        max = first(maximum(element_coords, dims=1))

        for p in vec_points
            # If point is outside the box check and the element check the next point
            all(p .>= min) && all(p .<= max) || continue
            # Check if the point is inside the element
            # p ⊄ e || continue

            # Main.@infiltrate
            weights = interpolate(e, p)
            # Main.@infiltrate
            for dof_symbol in dof_symbols
                dof_symbol_indexes = index.(dofs(e)[dof_symbol])
            end


            dofs, shape_functions = interpolate(e, p)

        end


        #         element_dofs = dofs(e)

        #         for (dof_symbol, dof_indexes) in pairs(element_dofs)
        #             Main.@infiltrate
        #             if haskey(interpolation_matrices, dof_symbol)

        #             end
        #         end
        #         # dofs, shape_functions = interpolate(e, p)


        #         # Main.@infiltrate

        #     end
        # end

    end

    return nothing
end

