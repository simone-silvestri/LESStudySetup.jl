using MPI
MPI.Init()
using LESStudySetup
using LESStudySetup.Oceananigans.Units
using JLD2

# Architecture (CPU, GPU, or Distributed)
arch = Distributed(GPU(), partition = Partition(4))

# Setting some initial values (Q = heat flux in W/m², Δz = vertical spacing)
set_value!(Δh = 4, Δz = 2, Lx = 10kilometers, Ly = 10kilometers)
set_value!(; Q = 50, τw = 0.0)

# Show all the parameters we are using
@info "Simulation parameters: " parameters

stop_time = 10days

# Let's start with an nonhydrostatic setup running for 30 days
simulation = idealized_setup(arch; stop_time)

if arch.local_rank == 0
    jldsave("experiment_$(experiment)_metadata.jld2", parameters = parameters)
end

# Show the configuration of the simulation
@info simulation

# Let's attach some outputs
model         = simulation.model
output_fields = merge(model.velocities, model.tracers)

simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                            schedule = TimeInterval(3hours),
                                                            overwrite_existing = true,
                                                            filename = "nonhydrostatic_snapshots_$(experiment)_$(arch.local_rank)")

#####
##### Let's run!!!!
#####

run!(simulation)