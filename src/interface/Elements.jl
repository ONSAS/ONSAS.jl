######################
# Elements interface #
######################

"""
Module defining the elements implemented.
"""
module Elements

@reexport import ..Utils: id

export Dof, id, symbol, is_fixed

""" Degree of freedom struct.
### Fields:
- `is_free` -- boolean indicating if the dof is free(`false`) or fixed(`true`).
- `symbol` -- degree of freedom symbol.
- `id`     -- degree of freedom id. 
"""
struct Dof
    is_fixed::Bool
    symbol::Symbol
    id::Int
end

"Returns the degree of freedom identification number."
id(d::Dof) = d.id

"Returns the degree of freedom symbol."
symbol(d::Dof) = d.symbol

"Returns `true` if the degree of freedom is free"
is_fixed(d::Dof) = d.is_fixed

"Set the degree of freedom as fixed"
fix!(d::Dof) = d.is_fixed = true

abstract type AbstractElement{M,D} end

""" Abstract supertype for all elements.

An `AbstractElement` object facilitates the process of evaluating:

    - The internal forces vector and its tangent matrices.
    - The inertial forces vector and its tangent matrices.
    - The mechanical stresses and loads.

    
**Common methods:**

* [`coordinates`](@ref)
* [`num_nodes`](@ref)
* [`dofs`](@ref)
* [`dofs_per_node`](@ref)
* [`loads`](@ref)
* [`loads_per_node`](@ref)
* [`geo_properties`](@ref)
* [`internal_force`](@ref)
* [`inertial_matrix`](@ref)
* [`mass_matrix`](@ref)
* [`nodes`](@ref)
* [`set_material!`](@ref)
* [`set_bc!`](@ref)
* [`set_nodes!`](@ref)
* [`solicitations`](@ref)

"""




end


