include("hydrostatic_experiment.jl")

run_experiment!("no_cooling_no_wind";   Q = 0.0,  τw = 0.0)
run_experiment!("cooling_10_wind_0075_σ_030"; Q = 10.0, τw = 0.075, σ² = 0.3)
run_experiment!("no_cooling_no_wind_no_restoring";          Q = 0.0,  τw = 0.0,   restoring = false)
