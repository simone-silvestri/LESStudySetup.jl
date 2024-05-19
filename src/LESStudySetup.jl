module LESStudySetup

export parameters
export idealized_setup
export set_value!, set!

using Reexport
@reexport using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Printf

include("xin_kai_vertical_diffusivity.jl")
include("parameters.jl")
include("utils.jl")
# include("initial_conditions.jl")
include("four_vortex_initial_conditions.jl")
include("idealized_setup.jl")
include("Diagnostics/Diagnostics.jl")

using .Diagnostics

end