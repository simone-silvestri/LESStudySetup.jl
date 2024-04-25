using LESStudySetup
using LESStudySetup.Oceananigans.Units

# Architecture (CPU, GPU, or Distributed)
architecture = CPU()

# Setting some initial values (Q = heat flux in W/m², Δz = vertical spacing)
set_value!(Q = 100, Δh = 200, Δz = 2)

# Show all the parameters we are using
@info "Simulation parameters: " parameters

# Let's start with an hydrostatic setup running for 30 days
stop_time  = 30days
simulation = idealized_setup(architecture; stop_time)

# Show the configuration of the simulation
@info simulation

# Let's attach some outputs
model         = simulation.model
output_fields = merge(model.velocities, model.tracers)

simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                         schedule = TimeInterval(3hours),
                                                         overwrite_existing = true,
                                                         filename = "nonhydrostatic_snapshots")

#####
##### Let's run!!!!
#####

run!(simulation)