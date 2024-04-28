using Oceananigans.Grids: node

model_type(::Val{true})  = HydrostaticFreeSurfaceModel
model_type(::Val{false}) = NonhydrostaticModel

# Forcing functions for the HydrostaticFreeSurfaceModel
@inline function Ur(i, j, k, grid, clock, fields, p) 
    x, y, z = node(i, j, k, grid, Face(), Center(), Center())
    return 1 / p.λ * (fields.u[i, j, k] - U̅(x, y, z, p))
end

@inline function Vr(i, j, k, grid, clock, fields, p) 
    x, y, z = node(i, j, k, grid, Center(), Face(), Center())
    return 1 / p.λ * (fields.v[i, j, k] - V̅(x, y, z, p))
end

@inline function Tr(i, j, k, grid, clock, fields, p) 
    x, y, z = node(i, j, k, grid, Center(), Face(), Center())
    return 1 / p.λ * (fields.T[i, j, k] - T̅(x, y, z, p))
end

function model_settings(::Type{HydrostaticFreeSurfaceModel}, grid)
    
    momentum_advection = WENO(; order = 7)
    tracer_advection   = WENO(; order = 7)
    closure = CATKEVerticalDiffusivity()
    tracers = (:T, :e)

    Fu = Forcing(Ur, discrete_form=true, parameters = merge((; λ = 5days), tuplify(parameters)))
    Fv = Forcing(Vr, discrete_form=true, parameters = merge((; λ = 5days), tuplify(parameters)))
    FT = Forcing(Tr, discrete_form=true, parameters = merge((; λ = 5days), tuplify(parameters)))

    Δh  = parameters.Δh 
    Lz  = parameters.Lz
    g   = parameters.g

    maximum_speed = 2 # m / s
    maximum_Δt    = 0.25 * Δh /  maximum_speed  

    wave_speed    = sqrt(Lz * g)
    baroclinic_Δt = 0.75 * Δh / wave_speed

    substeps = 2 * ceil(Int, maximum_Δt / baroclinic_Δt)

    @info "running with $substeps substeps"
    free_surface = SplitExplicitFreeSurface(grid; substeps, gravitational_acceleration = g)

    forcing = (; u = Fu, v = Fv, T = FT)

    return (; free_surface, momentum_advection, tracer_advection, tracers, closure, forcing)
end

# Background fields for the NonhydrostaticModel
@inline Tb(x, y, z, t, p) = T̅(x, y, z, p)
@inline Ub(x, y, z, t, p) = U̅(x, y, z, p)
@inline Vb(x, y, z, t, p) = V̅(x, y, z, p)

function model_settings(::Type{NonhydrostaticModel}, grid) 
    advection = WENO(; order = 7)
    tracers = :T

    T = BackgroundField(Tb; parameters = tuplify(parameters))
    U = BackgroundField(Ub; parameters = tuplify(parameters))
    V = BackgroundField(Vb; parameters = tuplify(parameters))

    background_fields = (; u = U, v = V, T)

    return (; advection, tracers, background_fields)
end

# Initial conditions for both models
@inline Ti(x, y, z) = Tᵢ(x, y, z, gpuify(parameters))
@inline Ui(x, y, z) =  U̅(x, y, z, gpuify(parameters))
@inline Vi(x, y, z) =  V̅(x, y, z, gpuify(parameters))

set_model!(model::HydrostaticFreeSurfaceModel) = 
    set!(model, u = Ui, v = Vi, T = Ti, e = 1e-6)

set_model!(model::NonhydrostaticModel) = set!(model, T = Ti)

function progress(sim) 
    u, v, w = sim.model.velocities
    T = sim.model.tracers.T

    ui = interior(u)
    vi = interior(v)
    wi = interior(w)

    Ti = interior(T)
    ei = interior(e)

    msg0 = @sprintf("Time: %s, iteration: %d, Δt: %s ", prettytime(sim.model.clock.time), 
                                                        sim.model.clock.iteration,
                                                        prettytime(sim.Δt))
    msg1 = @sprintf("(u, v, w): %.2e %.2e %.2e ", maximum(ui), maximum(vi), maximum(wi))
    msg2 = @sprintf("T: %.2e %.2e ", minimum(Ti), maximum(Ti))

    if sim.model isa HydrostaticFreeSurfaceModel
        e = sim.model.tracers.e
        msg3 = @sprintf("e: %.2e %.2e ", minimum(ei), maximum(ei))
    else
        msg3 = ""
    end

    @info msg0 * msg1 * msg2 * msg3 
end