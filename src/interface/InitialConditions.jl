#################################
# Initial Conditions interface #
#################################

"""
Module defining the boundary conditions implemented.
"""
module InitialConditions

export AbstractInitialCondition

""" Abstract supertype for all initial conditions type.

An `AbstractInitialCondition` object facilitates the process of defining:

    - Kinematic initial conditions.
    - Pre-strain initial conditions

**Common methods:**

* [`label`](@ref)
* [`set_label!`](@ref)
"""

abstract type AbstractInitialCondition end


end