module LESStudySetup

export parameters
export idealized_setup
export set_value!, set!

using Reexport
@reexport using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Printf

include("parameters.jl")
include("utils.jl")
include("initial_conditions.jl")
include("idealized_setup.jl")

end