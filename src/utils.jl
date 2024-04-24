using Oceananigans.Grids: node

model_type(::Val{true})  = HydrostaticFreeSurfaceModel
model_type(::Val{false}) = NonhydrostaticModel

function model_settings(::Type{HydrostaticFreeSurfaceModel}, grid)
    
    momentum_advection = WENO(; order = 7)
    tracer_advection   = WENO(; order = 7)
    closure = CATKEVerticalDiffusivity()
    tracers = (:T, :e)

    @inline function Ur(i, j, k, grid, clock, fields, p) 
        x, y, z = node(i, j, k, grid, Face(), Center(), Center())
        return 1 / p.λ * (fields.u[i, j, k] - U̅(x, y, z))
    end

    @inline function Vr(i, j, k, grid, clock, fields, p) 
        x, y, z = node(i, j, k, grid, Center(), Face(), Center())
        return 1 / p.λ * (fields.v[i, j, k] - V̅(x, y, z))
    end

    Fu = Forcing(Ur, discrete_form=true, parameters = (; λ = 10days))
    Fv = Forcing(Vr, discrete_form=true, parameters = (; λ = 10days))

    Δh  = parameters.Δh 
    Lz  = parameters.Lz
    g   = parameters.g

    maximum_speed = 20 # m / s
    maximum_Δt    = 0.25 * Δh /  maximum_speed  

    wave_speed    = sqrt(Lz * g)
    baroclinic_Δt = 0.75 * Δh / wave_speed

    substeps = 2 * ceil(Int, maximum_Δt / baroclinic_Δt)

    @info "running with $substeps substeps"
    free_surface = SplitExplicitFreeSurface(grid; substeps, gravitational_acceleration = g)

    forcing = (; u = Fu, v = Fv)

    return (; free_surface, momentum_advection, tracer_advection, tracers, forcing, closure)
end

function model_settings(::Type{NonhydrostaticModel}, grid) 
    advection = WENO(; order = 7)
    tracers = :T

    @inline Tb(x, y, z, t) = T̅(x, y, z)
    @inline Ub(x, y, z, t) = U̅(x, y, z)
    @inline Vb(x, y, z, t) = V̅(x, y, z)

    T = BackgroundField(Tb)
    U = BackgroundField(Tb)
    V = BackgroundField(Tb)

    background_fields = (; u = U, v = V, T)

    return (; advection, tracers, background_fields)
end

set_model!(model::HydrostaticFreeSurfaceModel) = set!(model, u = U̅, v = V̅, T = Tᵢ)
set_model!(model::NonhydrostaticModel) = set!(model, T = Tᵢ)

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