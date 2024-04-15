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
    ΔT :: Float64 = 1.0
    ρ₀ :: Float64 = 1020
    T₀ :: Float64 = 5
    cp :: Float64 = 3900
    N² :: Float64 = 2e-6
     H :: Float64 = 40
    ΔH :: Float64 = 15
    Δh :: Float64 = 1kilometers
    Δz :: Float64 = 4meters
    Lx :: Float64 = 100kilometers
    Ly :: Float64 = 200kilometers
    Lz :: Float64 = 250meters
     f :: Float64 = 1e-4
    τw :: Float64 = 0.1
     θ :: Float64 = 30
     Q :: Float64 = 0
     α :: Float64 = 2e-4
    Lf :: Float64 = 2
    σ² :: Float64 = 0.15
     g :: Float64 = Oceananigans.BuoyancyModels.g_Earth
end

Base.show(io::IO, c::ProblemConstants) =
    print(io, "├── temperature difference: ΔT = ", c.ΔT, "\n",
              "├── reference density:      ρ₀ = ", c.ρ₀, "\n",
              "├── surface temperature:    T₀ = ", c.T₀, "\n",
              "├── heat capacity:          cp = ", c.cp, "\n",
              "├── initial stratification  N² = ", c.N², "\n",
              "├── thermocline depth        H = ", c.H, "\n",
              "├── thermocline extent      ΔH = ", c.ΔH, "\n",
              "├── horizontal spacing:     Δh = ", c.Δh, "\n",
              "├── vertical spacing:       Δz = ", c.Δz, "\n",
              "├── x-domain size:          Lx = ", c.Lx, "\n",
              "├── y-domain size:          Ly = ", c.Ly, "\n",
              "├── z-domain size:          Lz = ", c.Lz, "\n",
              "├── Coriolis parameter:      f = ", c.f,  "\n",
              "├── wind stress:            τw = ", c.τw, "\n",
              "├── wind angle               θ = ", c.θ,  "\n",
              "├── heat flux:               Q = ", c.Q,  "\n", 
              "├── thermal expansion:       α = ", c.α,  "\n",
              "├── gravity:                 g = ", c.g,  "\n",
              "├── Initial vortex spread:  σ² = ", c.σ², "\n",
              "└── Frontal width:          Lf = ", c.Lf, "\n")

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
