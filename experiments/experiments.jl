include("hydrostatic_experiment.jl")

using Statistics: mean

run_experiment!("free_surface_short_test_050_wind_01"; Q = 50.0, τw = 0.1, Δh = 100)
run_experiment!("free_surface_short_test_050_wind_02"; Q = 50.0, τw = 0.2, Δh = 100)
run_experiment!("free_surface_short_test_000_wind_00"; Q = 00.0, τw = 0.0, Δh = 100)
