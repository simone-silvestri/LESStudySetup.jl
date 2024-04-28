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
                         N² = 5e-6)
    
    set_value!(; Q, τw, θ, ΔT, Lf, N²)

    @info "Simulation parameters: " parameters

    # Let's start with an hydrostatic setup running for 30 days
    stop_time  = 30days
    simulation = idealized_setup(architecture; stop_time, hydrostatic_approximation = true)

    jldsave("experiment_$(experiment)_metadata.jld2", parameters = parameters)

    # Show the configuration of the simulation
    @info simulation

    # Let's attach some outputs
    model         = simulation.model
    output_fields = merge(model.velocities, model.tracers)

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                             schedule = TimeInterval(3hours),
                                                             overwrite_existing = true,
                                                             filename = "hydrostatic_snapshots_exp$(experiment)")

    #####
    ##### Let's run!!!!
    #####

    run!(simulation)
end

# Quiescent experiment
run_experiment!("quiescent_ΔT2"; ΔT = 2)
run_experiment!("quiescent_ΔT4"; ΔT = 4)
run_experiment!("quiescent_ΔT6"; ΔT = 6)

# Weak Cooling Weak Wind
run_experiment!("weak_cooling_weak_wind_zonal"; Q = 50.0, τw = 0.075, θ = 0.0)
run_experiment!("weak_cooling_weak_wind_mixed"; Q = 50.0, τw = 0.075, θ = 45.0)
run_experiment!("weak_cooling_weak_wind_meridional"; Q = 50.0, τw = 0.075, θ = 90.0)

# Weak Cooling Strong Wind
run_experiment!("weak_cooling_strong_wind_zonal"; Q = 50.0, τw = 0.25, θ = 0.0)
run_experiment!("weak_cooling_strong_wind_mixed"; Q = 50.0, τw = 0.25, θ = 45.0)
run_experiment!("weak_cooling_strong_wind_meridional"; Q = 50.0, τw = 0.25, θ = 90.0)

# Strong Cooling Weak Wind
run_experiment!("strong_cooling_weak_wind_zonal"; Q = 350.0, τw = 0.075, θ = 0.0)
run_experiment!("strong_cooling_weak_wind_mixed"; Q = 350.0, τw = 0.075, θ = 45.0)
run_experiment!("strong_cooling_weak_wind_meridional"; Q = 350.0, τw = 0.075, θ = 90.0)

# Strong Cooling Strong Wind
run_experiment!("strong_cooling_strong_wind_zonal"; Q = 350.0, τw = 0.25, θ = 0.0)
run_experiment!("strong_cooling_strong_wind_mixed"; Q = 350.0, τw = 0.25, θ = 45.0)
run_experiment!("strong_cooling_strong_wind_meridional"; Q = 350.0, τw = 0.25, θ = 90.0)