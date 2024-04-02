module LESStudySetup

export parameters
export idealized_setup
export set_value!, set!

using Reexport
@reexport using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Printf

import Oceananigans.Fields: set!

using Base
using Base: getproperty

"""
    struct ProblemConstants

A mutable struct representing the constants used in the LES study setup.

# Fields
- `Δρ::Float64`: Density difference at the fronts.
- `ρ₀::Float64`: Reference density of the fluid.
- `N²::Float64`: Initial Buoyancy frequency squared.
- `Δh::Float64`: Horizontal grid spacing.
- `Δz::Float64`: Vertical grid spacing.
- `Lx::Float64`: Zonal domain size.
- `Ly::Float64`: Meridional domain size.
- `Lz::Float64`: Vertical domain depth.
- `Lf::Float64`: Width of the front.
- `θ::Float64`:  Vertical gradient of potential temperature.
- `f::Float64`:  Coriolis parameter.
- `τw::Float64`: Surface stress.
- `α::Float64`:  Thermal expansion coefficient.
- `Jᵀ::Float64`: Heat flux at the top.
"""
@kwdef mutable struct ProblemConstants
    Δρ :: Float64 = 0.25
    ρ₀ :: Float64 = 1029
    N² :: Float64 = 5e-6 
    Δh :: Float64 = 1kilometers
    Δz :: Float64 = 10meters
    Lx :: Float64 = 100kilometers
    Ly :: Float64 = 200kilometers
    Lz :: Float64 = 200meters
    Lf :: Float64 = 8kilometers
    cp :: Float64 = 3946
     θ :: Float64 = 0.25e-4
     f :: Float64 = 1e-4
    τw :: Float64 = 0.1
     α :: Float64 = 2e-4
     Q :: Float64 = 0
end

Base.show(io::IO, c::ProblemConstants) =
    print(io, "├── density difference:     Δρ = ", c.Δρ, "\n",
              "├── reference density:      ρ₀ = ", c.ρ₀, "\n",
              "├── heat capacity:          cp = ", c.cp, "\n",
              "├── initial stratification: N² = ", c.N², "\n",
              "├── horizontal spacing:     Δh = ", c.Δh, "\n",
              "├── vertical spacing:       Δz = ", c.Δz, "\n",
              "├── x-domain size:          Lx = ", c.Lx, "\n",
              "├── y-domain size:          Ly = ", c.Ly, "\n",
              "├── z-domain size:          Lz = ", c.Lz, "\n",
              "├── frontal width:          Lf = ", c.Lf, "\n",
              "├── maximum confluence       θ = ", c.θ,  "\n",
              "├── Coriolis parameter:      f = ", c.f,  "\n",
              "├── wind stress:            τw = ", c.τw, "\n",
              "├── thermal expansion:       α = ", c.α,  "\n",
              "└── heat flux:               Q = ", c.Q, "\n")

# The constants of the idealized setup
const parameters = ProblemConstants()

# Setting problem constants
set_value!(; kwargs...)      = set!(parameters; kwargs...)
set_value!(var::Symbol, val) = parameters[var] = val

Base.getindex(c::ProblemConstants, var::Symbol)     = @eval $c.$var
Base.setindex!(c::ProblemConstants, v, var::Symbol) = @eval $c.$var = $v

function set!(c::ProblemConstants; kwargs...)
    for (fldname, value) in kwargs
        if fldname ∈ propertynames(c)
            c[fldname] = value
        end
    end

    return nothing
end

include("utils.jl")
include("initial_conditions.jl")
include("idealized_setup.jl")

end