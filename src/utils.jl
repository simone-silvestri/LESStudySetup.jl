using Oceananigans.Grids: node
using Statistics: mean

model_type(::Val{true})  = HydrostaticFreeSurfaceModel
model_type(::Val{false}) = NonhydrostaticModel

# Forcing functions for the HydrostaticFreeSurfaceModel
@inline function Ur(i, j, k, grid, clock, fields, p) 
    x, y, z = node(i, j, k, grid, Face(), Center(), Center())
    return @inbounds 1 / p.λ * (U̅(x, y, z, p) - p.Um[i, j, k])
end

@inline function Vr(i, j, k, grid, clock, fields, p) 
    x, y, z = node(i, j, k, grid, Center(), Face(), Center())
    return @inbounds 1 / p.λ * (V̅(x, y, z, p) - fields.v[i, j, k])
end

@inline function Tr(i, j, k, grid, clock, fields, p) 
    x, y, z = node(i, j, k, grid, Center(), Center(), Center())
    return @inbounds 1 / p.λ * (T̅(x, y, z, p) - p.Tm[i, j, k])
end

function model_settings(::Type{HydrostaticFreeSurfaceModel}, grid;
                        restoring = false)
    
    momentum_advection = WENO(; order = 7)
    tracer_advection   = WENO(; order = 7)
    closure = XinKaiVerticalDiffusivity()
    tracers = :T

    Tm = Field{Center, Nothing, Nothing}(grid)
    Um = Field{Face,   Center,  Nothing}(grid)
    
    Fu = Forcing(Ur, discrete_form=true, parameters = merge((; λ = 5days), tuplify(parameters), (; Um)))
    Fv = Forcing(Vr, discrete_form=true, parameters = merge((; λ = 5days), tuplify(parameters)))
    FT = Forcing(Tr, discrete_form=true, parameters = merge((; λ = 5days), tuplify(parameters), (; Tm)))

    free_surface = SplitExplicitFreeSurface(grid; substeps = 75, gravitational_acceleration = parameters.g)
    @info "running with $(length(free_surface.settings.substepping.averaging_weights)) substeps"

    forcing = if restoring
        (; u = Fu, v = Fv, T = FT)
    else
        NamedTuple()
    end

    return (; free_surface, momentum_advection, tracer_advection, tracers, closure, forcing)
end

# Background fields for the NonhydrostaticModel
@inline Tb(x, y, z, t) = Tᵢ(x, y, z)
@inline Ub(x, y, z, t) = uᵢ(x, y, z)
@inline Vb(x, y, z, t) = vᵢ(x, y, z)

function model_settings(::Type{NonhydrostaticModel}, grid;
                        restoring = false) 
    advection = WENO(; order = 7)
    tracers = :T

    T = BackgroundField(Tb; parameters = tuplify(parameters))
    U = BackgroundField(Ub; parameters = tuplify(parameters))
    V = BackgroundField(Vb; parameters = tuplify(parameters))

    background_fields = if restoring
        (; u = U, v = V, T)
    else
        NamedTuple()
    end

    return (; advection, tracers, background_fields)
end

# Initial conditions for both models
@inline Ti(x, y, z) = Tᵢ(x, y, z)
@inline Ui(x, y, z) = uᵢ(x, y, z)
@inline Vi(x, y, z) = vᵢ(x, y, z)

set_model!(model::HydrostaticFreeSurfaceModel) = 
    set!(model, u = uᵢ, v = vᵢ, T = Tᵢ)

set_model!(model::NonhydrostaticModel) = set!(model, T = Tᵢ)

simulation_callbacks!(::Simulation{<:NonhydrostaticModel}, args...) = nothing

function simulation_callbacks!(sim::Simulation{<:HydrostaticFreeSurfaceModel}, restoring) 
    if !restoring
        return nothing
    end

    function update_mean!(sim)
        sim.model.forcing.T.parameters.Tm .= mean(sim.model.tracers.T,    dims = (2, 3))
        sim.model.forcing.u.parameters.Um .= mean(sim.model.velocities.u, dims = 3)
    end

    sim.callbacks[:update_mean] = Callback(update_mean!, IterationInterval(10))

    return nothing
end

function progress(sim) 
    u, v, w = sim.model.velocities
    T = sim.model.tracers.T

    ui = interior(u)
    vi = interior(v)
    wi = interior(w)

    Ti = interior(T)

    msg0 = @sprintf("Time: %s, iteration: %d, Δt: %s ", prettytime(sim.model.clock.time), 
                                                        sim.model.clock.iteration,
                                                        prettytime(sim.Δt))
    msg1 = @sprintf("(u, v, w): %.2e %.2e %.2e ", maximum(ui), maximum(vi), maximum(wi))
    msg2 = @sprintf("T: %.2e %.2e ", minimum(Ti), maximum(Ti))

    @info msg0 * msg1 * msg2 
end
