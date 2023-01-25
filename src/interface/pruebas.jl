abstract type AbstractLoads end

mutable struct NodalLoadsBoundaryCondition <: AbstractLoads
    loadsBaseVals::Vector
    loadsCoordSystem::String
    loadsTimeFactor::Float64
    function NodalLoadsBoundaryCondition(loadsBaseVals=nothing, loadsCoordSystem=nothing, loadsTimeFactor=nothing)
        new([], "", 1.0)
    end
end
# constructor with missing fields
function NodalLoadsBoundaryCondition(loadsBaseVals, loadsCoordSystem)
    return NodalLoadsBoundaryCondition(loadsBaseVals, loadsCoordSystem, 1.0)
end
# function NodalLoadsBoundaryCondition() end

struct UserLoadsBoundaryCondition <: AbstractLoads
    user_load_function::Function
    function UserLoadsBoundaryCondition(user_load_function=nothing)
        function f() end
        new(f)
    end
end

struct ElementLoads <: AbstractLoads ## podria ser uniform element loads, otro struct para variables, etc
    # to do
    q
    function ElementLoads(q=nothing)
        new(q)
    end
end

struct LoadsBoundaryCondition <: AbstractLoads
    nodal_loads::NodalLoadsBoundaryCondition
    user_loads::UserLoadsBoundaryCondition
    element_loads::ElementLoads
end

function LoadsBoundaryCondition(nodal::NodalLoadsBoundaryCondition)
    function f() end
    return LoadsBoundaryCondition(nodal, UserLoadsBoundaryCondition(), ElementLoads())
end
function LoadsBoundaryCondition(nodal::NodalLoadsBoundaryCondition, user::UserLoadsBoundaryCondition)
    return LoadsBoundaryCondition(nodal, user, ElementLoads())
end
function LoadsBoundaryCondition(nodal::NodalLoadsBoundaryCondition, element::ElementLoads)
    return LoadsBoundaryCondition(nodal, UserLoadsBoundaryCondition(), element)
end

a = NodalLoadsBoundaryCondition(ones(6), "global")
b = UserLoadsBoundaryCondition()
c = ElementLoads()

d = LoadsBoundaryCondition(a, b, c)