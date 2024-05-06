using LESStudySetup
using LESStudySetup.Oceananigans.Units
using JLD2 

# Architecture (CPU, GPU, or Distributed)
architecture = GPU()

# Setting some initial values (Q = heat flux in W/m², Δz = vertical spacing)
set_value!(Δh = 200, Δz = 2)

function run_experiment!(experiment; 
                         Q  = 0.0, # Cooling heat flux in W/m²
                         τw = 0.0, # Wind stress in N/m²
                         θ  = 0.0, # Wind stress angle in degrees (0 correspond to zonal wind stress)
                         ΔT = 2.0, # Meridional temperature difference
                         Lf = 1.0, # Size of temperature front (large numbers correspond to steeper fronts)
                         N² = 5e-6,
                         restoring = true)
    
    set_value!(; Q, τw, θ, ΔT, Lf, N²)

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
                                                             schedule = TimeInterval(3hours),
                                                             overwrite_existing = true,
                                                             array_type = Array{Float32},
                                                             filename = "hydrostatic_snapshots_$(experiment)")

    #####
    ##### Let's run!!!!
    #####

    run!(simulation)

    return simulation
end

run_experiment!("constant_cooling_05_wind_0075_30"; Q = 5.0,  τw = 0.075, θ = 30.0)
run_experiment!("constant_cooling_10_wind_0075_30"; Q = 10.0, τw = 0.075, θ = 30.0)
run_experiment!("constant_cooling_25_wind_0075_30"; Q = 25.0, τw = 0.075, θ = 30.0)
run_experiment!("constant_cooling_50_wind_0075_30"; Q = 50.0, τw = 0.075, θ = 30.0)
run_experiment!("no_cooling_wind_0075_30";          Q = 0.0,  τw = 0.075, θ = 30.0)
run_experiment!("no_cooling_no_wind";               Q = 0.0,  τw = 0.0,   θ = 30.0)

# Experiments without restoring, only wind and cooling
run_experiment!("constant_cooling_10_wind_0075_30_no_restoring"; Q = 10.0, τw = 0.075, θ = 30.0, restoring = false)




