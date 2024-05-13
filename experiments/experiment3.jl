include("hydrostatic_experiment.jl")

run_experiment!("cooling_10_wind_0075_no_restoring";        Q = 10.0, τw = 0.075, restoring = false)
run_experiment!("cooling_10_wind_0075_σ_030_no_restoring";  Q = 10.0, τw = 0.075, σ² = 0.3, restoring = false)


