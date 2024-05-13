using LESStudySetup
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Oceananigans.Utils: ConsecutiveIterations
using JLD2 

# Architecture (CPU, GPU, or Distributed)
architecture = GPU()

# Setting some initial values (Q = heat flux in W/m², Δz = vertical spacing)
set_value!(Δh = 250, Δz = 2)

function run_experiment!(experiment; 
                         Q  = 0.0,  # Cooling heat flux in W/m²
                         τw = 0.0,  # Wind stress in N/m²
                         θ  = 30.0, # Wind stress angle in degrees (0 correspond to zonal wind stress)
                         ΔT = 2.0,  # Meridional temperature difference
                         Lf = 1.0,  # Size of temperature front (large numbers correspond to steeper fronts)
                         N² = 5e-6,
                         σ² = 0.15,
                         restoring = true)
    
    set_value!(; Q, τw, θ, ΔT, Lf, N², σ²)

    @info "Simulation parameters: " parameters

    # Let's start with an hydrostatic setup running for 20 days
    stop_time  = 20days
    simulation = idealized_setup(architecture; stop_time, hydrostatic_approximation = true, restoring)

    jldsave("experiment_$(experiment)_metadata.jld2", parameters = parameters)

    # Show the configuration of the simulation
    @info simulation

    # Let's attach some outputs
    model         = simulation.model
    output_fields = merge(model.velocities, model.tracers)

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                             schedule = ConsecutiveIterations(TimeInterval(3hours)),
                                                             overwrite_existing = true,
                                                             array_type = Array{Float32},
                                                             filename = "hydrostatic_snapshots_$(experiment)")

    #####
    ##### Let's run!!!!
    #####

    run!(simulation)

    return simulation
end

# run_experiment!("no_cooling_wind_0075"; Q = 0.0,  τw = 0.075)
# run_experiment!("cooling_05_wind_0075"; Q = 5.0,  τw = 0.075)
# run_experiment!("cooling_10_wind_0075"; Q = 10.0, τw = 0.075)
# run_experiment!("cooling_35_wind_0075"; Q = 35.0, τw = 0.075)
# run_experiment!("no_cooling_no_wind";   Q = 0.0,  τw = 0.0)
# 
# 
# run_experiment!("cooling_10_wind_0075_σ_030"; Q = 10.0, τw = 0.075, σ² = 0.3)
# 
# # Experiments without restoring
# run_experiment!("no_cooling_no_wind_no_restoring";          Q = 0.0,  τw = 0.0,   restoring = false)
# run_experiment!("cooling_10_wind_0075_no_restoring";        Q = 10.0, τw = 0.075, restoring = false)
# run_experiment!("cooling_10_wind_0075_σ_030_no_restoring";  Q = 10.0, τw = 0.075, σ² = 0.3, restoring = false)




