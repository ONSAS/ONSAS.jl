"""
Module defining structural analyses that can be solved. 
"""
module StructuralAnalyses

using ..StructuralModel

export StaticAnalysis, structure, final_time, time_step


""" Abstract supertype for all structural analysis.

An `AbstractStructuralAnalysis` object facilitates the process of defining an structural analysis
to be solved.

**Common methods:**

* [`structure`](@ref)
* [`final_time`](@ref)
* [`time_step`](@ref)
"""
abstract type AbstractStructuralAnalysis end

"Returns analyzed structure"
structure(a::AbstractStructuralAnalysis) = a.s

"Returns final time of the analysis"
final_time(a::AbstractStructuralAnalysis) = a.t₁

"Returns time step of the analysis"
time_step(a::AbstractStructuralAnalysis) = a.Δₜ

# ================
# Static analysis
# ================

""" StaticAnalysis struct.
A `StaticAnalysis` is a collection of parameters for defining the static analysis of the structure. 

### Fields:
- `s` -- Stores the structure to be analyzed.
- `t₁` -- Stores the final time of the analysis.
- `t₁` -- Stores the final time of the analysis.
- `initial_state` -- Stores the initial state of the structure.
"""

struct StaticAnalysis <: AbstractStructuralAnalysis
    s::AbstractStructure
    t₁::Number
    Δₜ::Number
    init_state::StaticState
    function StaticAnalysis(s, t₁::Number=1.0, Δₜ::Number=0.1, init_state=current_state(s))
        # Check meaningful parameters
        (t₁ > 0 && Δₜ > 0) || throw(ArgumentError("t₁ and Δₜ must be positive"))
        (t₁ ≥ Δₜ) || throw(ArgumentError("t₁ must be greater than Δₜ"))
        new(s, t₁, Δₜ, init_state)
    end
end

# ================
# Dynamic analysis
# ================



# ================
# Modal analysis
# ================



end #module 