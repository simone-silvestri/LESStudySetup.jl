using LESStudySetup
using LESStudySetup.Oceananigans.Units

# Architecture (CPU, GPU, or Distributed)
architecture = GPU()

# Setting some initial values (Q = heat flux in W/m², Δz = vertical spacing)
set_value!(Q = 100, Δh = 20000, Δz = 20)

# all_Q  = [0, 50, 100]
# all_τw = [0, 0.1, 0.2]
# all_θ  = [0, 45, 90]
# all_ΔT = [1, 3, 6]
# all_Lf = [1, 2, 3]
# all_N² = [1e-6, 5e-6, 1e-5]

# experiment = 1

# for Q in all_Q, τw in all_τw, θ in all_θ, ΔT in all_ΔT, Lf in all_Lf, N² in all_N²

#     set_value!(; Q, τw, θ, ΔT, Lf, N²)

    # Show all the parameters we are using
    @info "Simulation parameters: " parameters

    # Let's start with an hydrostatic setup running for 30 days
    stop_time  = 30days
    simulation = idealized_setup(architecture; stop_time, hydrostatic_approximation = true)

    # Show the configuration of the simulation
    @info simulation

    # Let's attach some outputs
    model         = simulation.model
    output_fields = merge(model.velocities, model.tracers)

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                             schedule = TimeInterval(3hours),
                                                             overwrite_existing = true,
                                                             filename = "hydrostatic_snapshots")

    #####
    ##### Let's run!!!!
    #####

    run!(simulation)
# end