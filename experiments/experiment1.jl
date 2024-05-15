include("hydrostatic_experiment.jl")

using CUDA
CUDA.device!(1)
using Statistics: mean

# run_experiment!("no_cooling_wind_0075"; Q = 0.0,  τw = 0.075)
# run_experiment!("cooling_05_wind_0075"; Q = 5.0,  τw = 0.075)
# run_experiment!("cooling_10_wind_0075"; Q = 10.0, τw = 0.075)
run_experiment!("strong_stratification_cooling_100_wind_0075"; Q = 100.0, τw = 0.075, N² = 1.2e-5)
run_experiment!("strong_stratification_cooling_50_wind_0075";  Q = 35.0,  τw = 0.075, N² = 1.2e-5)
