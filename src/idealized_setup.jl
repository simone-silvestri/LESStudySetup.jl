model_type(::Val{true})  = HydrostaticFreeSurfaceModel
model_type(::Val{false}) = NonhydrostaticModel

model_specific_kwargs(::Val{true})  = (; momentum_advection = WENO(; order = 7),
                                         tracer_advection = WENO(; order = 7),
                                         closure = CATKEVerticalDiffusivity(),
                                         tracers = (:T, :e))

model_advection(::Val{false}) = (; advection = WENO(; order = 7), 
                                   tracers = :T)

function progress(sim) 
    u, v, w = sim.model.velocities
    T = sim.model.tracers.T

    msg0 = @sprintf("Time: %s, iteration: %d, Δt: %s ", prettytime(sim.model.clock.time), 
                                                        sim.model.clock.iteration,
                                                        prettytime(sim.Δt))
    msg1 = @sprintf("(u, v, w): %.2e %.2e %.2e ", maximum(abs, u), maximum(abs, v), maximum(abs, w))
    msg2 = @sprintf("T: %.2e %.2e ", minimum(T), maximum(T))

    if sim.model isa HydrostaticFreeSurfaceModel
        e = sim.model.tracers.e
        msg3 = @sprintf("e: %.2e %.2e ", minimum(e), maximum(e))
    else
        msg3 = ""
    end

    @info msg0 * msg1 * msg2 * msg3 
end
"""
    idealized_setup(arch; stop_time = 100days, hydrostatic_approximation = true)

Create and configure a simulation for an idealized LES setup.

## Arguments
- `arch`: The architecture to use for the simulation.
- `stop_time`: The duration of the simulation in days. Default is 100 days.
- `hydrostatic_approximation`: Whether to use the hydrostatic approximation. Default is `true`.

## Returns
A `Simulation` object configured for the idealized setup.

The function retrieves the problem constants and calculates the grid size based on the problem dimensions. 
It then creates a `RectilinearGrid` object with the specified architecture and grid parameters. 
The model type and specific keyword arguments are determined based on the hydrostatic approximation setting. 
The function sets up the boundary conditions for velocity and temperature fields, 
and creates a `Model` object with the specified grid, coriolis, buoyancy, boundary conditions, and keyword arguments.
The initial conditions for velocity and temperature fields are set, and the maximum velocity magnitude is calculated to determine the time step size. 
A `TimeStepWizard` object is created to control the time step size during the simulation. 
Finally, a `Simulation` object is created with the model, time step, and stop time, and the progress and wizard callbacks are added.
"""
function idealized_setup(arch; 
                         stop_time = 100days,
                         hydrostatic_approximation = true)
    
    # Retrieving the problem constants
    Δh = problem_constants.Δh 
    Δz = problem_constants.Δz 
    Lx = problem_constants.Lx 
    Ly = problem_constants.Ly 
    Lz = problem_constants.Lz 
     α = problem_constants.α
     f = problem_constants.f  
    ρ₀ = problem_constants.ρ₀
    τw = problem_constants.τw 
    Jᵀ = problem_constants.Jᵀ

    Nx = ceil(Int, Lx / Δh)
    Ny = ceil(Int, Ly / Δh)
    Nz = ceil(Int, Lz / Δz)

    grid = RectilinearGrid(arch, 
                           size = (Nx, Ny, Nz), 
                           x = (0, Lx), 
                           y = (0, Ly), 
                           z = (0, Lz),
                           halo = (4, 4, 4))

    ModelType = model_type(Val(hydrostatic_approximation))
    kwargs    = model_specific_kwargs(Val(hydrostatic_approximation))

    coriolis = FPlane(; f)
    buoyancy = SeawaterBuoyancy(; equation_of_state = LinearEquationOfState(thermal_expansion = α), 
                                  constant_salinity = 35)
    
    u_top = FluxBoundaryCondition(τw / ρ₀)
    T_top = FluxBoundaryCondition(Jᵀ)

    u_bcs = FieldBoundaryConditions(top = u_top)
    T_bcs = FieldBoundaryConditions(top = T_top)

    boundary_conditions = (u = u_bcs, T = T_bcs)
    
    model = ModelType(; grid, 
                        coriolis,
                        buoyancy,
                        boundary_conditions,
                        kwargs...)

    set!(model, u = uᵢ, v = vᵢ, T = Tᵢ)

    u, v, w = model.velocities

    u_max = max(maximum(abs, u), maximum(abs, v))

    Δt = 0.2 * Δh / u_max
    
    wizard = TimeStepWizard(cfl = 0.25, max_change = 1.1)

    simulation = Simulation(model; Δt, stop_time)

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    return simulation
end
