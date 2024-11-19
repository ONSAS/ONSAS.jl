"""
Module to handle Dirichlet boundary conditions.
"""
module DirichletBoundaryConditions

using ..Utils
using ..BoundaryConditions

export Dirichlet

"""
Generalized displacement boundary condition.
"""
struct Dirichlet <: AbstractDirichletBoundaryCondition
    "Field where the boundary condition applies to."
    field::Field
    "Function defining values which are imposed by this boundary condition."
    values::Function
    "Boundary condition label."
    name::Label
    function Dirichlet(field::Field, values::Function, name::Label = NO_LABEL)
        new(field, values, name)
    end
end

end
