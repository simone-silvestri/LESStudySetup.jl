module LESStudySetup

using Oceananigans
using Oceananigans.Units

import Oceananigans.Fields: set!

@kwdef mutable struct ProblemConstants
   Δbʸ :: Float64 = 1e-8 
    ρ₀ :: Float64 = 1029
    N² :: Float64 = 5e-6 
    Δh :: Float64 = 1kilometer
    Lx :: Float64 = 100kilometers
    Ly :: Float64 = 100kilometers
    Lz :: Float64 = 200
    Lf :: Float64 = 100
     α :: Float64 = 0.25e-4
     f :: Float64 = 1e-4
    τw :: Float64 = 0.1
end

# The constants in the idealized setup
problem_constants = ProblemConstants()

# Setting problem constants
set_value!(; kwargs...) = set!(problem_constants; kwargs...)

# Need to change this
function set!(c::ProblemConstants; kwargs...)
    for (fldname, value) in kwargs
        if fldname ∈ propertynames(c)
            ϕ = getproperty(c, fldname)
            ϕ = value
        end
    end

    return nothing
end

end