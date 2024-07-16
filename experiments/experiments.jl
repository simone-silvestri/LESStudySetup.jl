include("hydrostatic_experiment.jl")

using Statistics: mean

run_experiment!("free_surface_short_test_100_wind_00"; Q = 50.0, τw = 0.1, Δh = 100)