"Module defining the interface to handle with different interpolators."
module Interpolators

using Dictionaries, Reexport

using ..Elements

export AbstractInterpolator, FEMInterpolator, points, interpolate, node_to_weights,
       points_to_element

"""
An `AbstractInterpolator ` is an object that stores the weights nodes and elements needed to interpolate
a magnitude at a given points.

*Common methods:*
- [`interpolate`](@ref)
"""
abstract type AbstractInterpolator end

"""
A `FEMInterpolator` struct stores the weights nodes and elements needed
to interpolate the solution at a given set of points. The weight are given 
by the element's shape functions. The index of each `Vector` is the index 
in the `Vector` of `Point`s where the magnitude is evaluated.
"""
struct FEMInterpolator{dim,
                       P<:Point{dim},
                       VP<:AbstractVector{P},
                       N<:AbstractNode{dim},
                       TW<:Real,
                       E<:AbstractElement,
                       VE<:AbstractVector{E}} <: AbstractInterpolator
    "`Vector` of `Point`s where the solution is interpolated."
    points_interpolated::VP
    "`Dictionary` with `Node`s as keys and the corresponding weights as values."
    node_to_weights::Vector{Dictionary{N,TW}}
    "`Element` where the point is located."
    points_to_element::VE
    # Check the lengths
    function FEMInterpolator(points_interpolated::VP,
                             node_to_weights::Vector{Dictionary{N,TW}},
                             points_to_element::VE) where {dim,
                                                           P<:Point{dim},
                                                           VP<:AbstractVector{P},
                                                           N<:AbstractNode,
                                                           TW<:Real,
                                                           E<:AbstractElement,
                                                           VE<:AbstractVector{E}}
        @assert length(points_interpolated) == length(node_to_weights) == length(points_to_element) "All FEMInterpolator inputs must have the same length"
        new{dim,P,VP,N,TW,E,VE}(points_interpolated, node_to_weights, points_to_element)
    end
end

"Return a `Vector` of `Point`s where the solution is evaluated."
points(fem_interpolator::FEMInterpolator) = fem_interpolator.points_interpolated

"Return a `Vector` that maps the weight corresponding to each `Node`."
node_to_weights(fem_interpolator::FEMInterpolator) = fem_interpolator.node_to_weights

"Return a `Vector` that maps the `Element` where each `Point` is located."
points_to_element(fem_interpolator::FEMInterpolator) = fem_interpolator.points_to_element

"Return the `Real` `nodal_magnitude` interpolated with a `FEMInterpolator`."
function interpolate(nodal_magnitude::Dictionary{<:AbstractNode{dim},<:Real},
                     fem_interpolator::FEMInterpolator{dim}) where {dim}
    # Weights for each node
    node_to_w = node_to_weights(fem_interpolator)
    num_points_to_interpolate = length(points(fem_interpolator))
    # Initialize the interpolated magnitude
    interpolated_magnitude = zeros(num_points_to_interpolate)
    # Loop over the points and sum up the contributions
    for point_idx in 1:num_points_to_interpolate
        for (node, weight) in pairs(node_to_w[point_idx])
            interpolated_magnitude[point_idx] += weight * nodal_magnitude[node]
        end
    end
    interpolated_magnitude
end

"Return the `Vector{Real}` `nodal_magnitude` interpolated with a `FEMInterpolator`."
function interpolate(nodal_magnitude::Dictionary{<:AbstractNode{dim},<:AbstractVector{<:Real}},
                     fem_interpolator::FEMInterpolator{dim}) where {dim}
    # Weights for each node
    node_to_w = node_to_weights(fem_interpolator)
    num_points_to_interpolate = length(points(fem_interpolator))
    # Dimension of the nodal magnitude vectors
    nodal_mag_dim = length(first(values(nodal_magnitude)))
    # Initialize the interpolated magnitude
    interpolated_magnitude = [zeros(nodal_mag_dim) for _ in 1:num_points_to_interpolate]
    # Loop over the points and sum up the contributions
    for point_idx in 1:num_points_to_interpolate
        for (node, weight) in pairs(node_to_w[point_idx])
            nodal_mag_components = nodal_magnitude[node]
            for i in 1:nodal_mag_dim
                interpolated_magnitude[point_idx][i] += weight * nodal_mag_components[i]
            end
        end
    end
    interpolated_magnitude
end

end
