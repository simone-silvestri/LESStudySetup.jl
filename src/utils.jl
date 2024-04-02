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