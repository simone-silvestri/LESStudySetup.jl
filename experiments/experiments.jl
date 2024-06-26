include("hydrostatic_experiment.jl")

using Statistics: mean

# run_experiment!("four_vortices_cooling_100_wind_00"; Q = 100.0, τw = 0.0)
# run_experiment!("four_vortices_cooling_010_wind_00"; Q = 10.0,  τw = 0.0)
# run_experiment!("four_vortices_cooling_050_wind_00"; Q = 50.0,  τw = 0.0)
# run_experiment!("four_vortices_cooling_100_wind_01"; Q = 100.0, τw = 0.1)
# run_experiment!("four_vortices_cooling_010_wind_01"; Q = 10.0,  τw = 0.1)
# run_experiment!("four_vortices_cooling_050_wind_01"; Q = 50.0,  τw = 0.1)

# for prefix in ["four_vortices_cooling_100_wind_00",
#                "four_vortices_cooling_010_wind_00",
#                "four_vortices_cooling_050_wind_00",
#                "four_vortices_cooling_100_wind_01",
#                "four_vortices_cooling_010_wind_01",
#                "four_vortices_cooling_050_wind_01"]

#     _ = write_pointwise_diagnostics(prefix; architecture = GPU())
# end

run_experiment!("free_surface_short_test_100_wind_00"; Q = 50.0, τw = 0.1, 
                                                       stop_time = 5hours,
                                                       output_frequency = 5minutes)