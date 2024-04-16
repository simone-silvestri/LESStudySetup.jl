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
    Δh = parameters.Δh 
    Δz = parameters.Δz 
    Lx = parameters.Lx 
    Ly = parameters.Ly 
    Lz = parameters.Lz 
     α = parameters.α
     f = parameters.f  
    ρ₀ = parameters.ρ₀
    cₚ = parameters.cp
    τw = parameters.τw 
     θ = parameters.θ
     Q = parameters.Q

    # Calculating the grid-size
    Nx = ceil(Int, Lx / Δh)
    Ny = ceil(Int, Ly / Δh)
    Nz = ceil(Int, Lz / Δz)

    # Constructing the grid
    grid = RectilinearGrid(arch, 
                           size = (Nx, Ny, Nz), 
                           x = (0, Lx), 
                           y = (0, Ly), 
                           z = (-Lz, 0),
                           halo = (4, 4, 4))

    # ModelType can be either a `HydrostaticFreeSurfaceModel` or a `NonhydrostaticModel`
    ModelType = model_type(Val(hydrostatic_approximation))
    kwargs    = model_specific_kwargs(Val(hydrostatic_approximation))

    coriolis = FPlane(; f)
    buoyancy = SeawaterBuoyancy(; equation_of_state = LinearEquationOfState(thermal_expansion = α), 
                                  constant_salinity = 35)
    
    # Cooling in the middle of the domain and heating outside?
    @inline Qtop(x, y, t, p) = - p.Q / p.ρ₀ / p.cₚ * cos(2π * x / p.Lx)

    u_top = FluxBoundaryCondition(τw * cosd(θ) / ρ₀)
    v_top = FluxBoundaryCondition(τw * sind(θ) / ρ₀)
    T_top = FluxBoundaryCondition(Qtop, parameters = (; Q, Lx, cₚ, ρ₀)) # Positive fluxes at the top are cooling in Oceananigans

    u_bcs = FieldBoundaryConditions(top = u_top)
    v_bcs = FieldBoundaryConditions(top = v_top)
    T_bcs = FieldBoundaryConditions(top = T_top)

    boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs)
    
    model = ModelType(; grid, 
                        coriolis,
                        buoyancy,
                        boundary_conditions,
                        kwargs...)

    set!(model, u = U̅, v = V̅, T = Tᵢ)

    u, v, w = model.velocities

    u_max = max(maximum(abs, u), maximum(abs, v))

    Δt = 0.2 * Δh / u_max
    
    wizard = TimeStepWizard(cfl = 0.25, max_change = 1.1)

    simulation = Simulation(model; Δt, stop_time)

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    return simulation
end
