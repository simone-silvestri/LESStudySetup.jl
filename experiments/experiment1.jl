include("hydrostatic_experiment.jl")

using Statistics: mean

run_experiment!("no_cooling_wind_0075"; Q = 0.0,  τw = 0.075)
run_experiment!("cooling_05_wind_0075"; Q = 5.0,  τw = 0.075)
run_experiment!("cooling_10_wind_0075"; Q = 10.0, τw = 0.075)
run_experiment!("cooling_35_wind_0075"; Q = 35.0, τw = 0.075)
