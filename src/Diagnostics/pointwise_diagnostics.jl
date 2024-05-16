using Oceananigans.Operators: div_xyᶜᶜᶜ

""" propagate a diagnostic over the timeseries found in snapshots 
    and save the results in a FieldTimeSeries object that is backed up in 
    `filename`.
"""
function propagate_function(func, snapshots; filename = "temp.jld2")
    first_operation = func(snapshots, 1)
    field = Field(first_operation)
    compute!(field)

    loc  = location(field)
    grid = field.grid

    saved_times = times(snapshots)
    func_name   = String(Symbol(func))

    field_time_series = FieldTimeSeries{loc...}(grid, Nt; 
                                                times = saved_times,
                                                backend = OnDisk(),
                                                path = filename,
                                                name = func_name)

    @info "calculating $func_name over the timeseries..."

    for i in 1:Nt
        set!(field, func(snapshots, i))
        set!(field_time_series, field, i)
        @info "calculated $func_name at iteration $i of $Nt"
    end

    return field_time_series
end

""" x-z momentum flux """
function uw(snapshots, i)
    u = snapshots[:u][i]
    w = snapshots[:w][i]

    return u * w
end

""" y-z momentum flux """
function vw(snapshots, i)
    v = snapshots[:u][i]
    w = snapshots[:w][i]

    return v * w
end

""" zonal buoyancy flux """
function ub(snapshots, i)
    α = parameters.α
    g = parameters.g

    u = snapshots[:u][i]
    T = snapshots[:T][i]
    
    return α * g * T * u
end

""" meridional buoyancy flux """
function vb(snapshots, i)
    α = parameters.α
    g = parameters.g

    v = snapshots[:v][i]
    T = snapshots[:T][i]
    
    return α * g * T * v
end

""" vertical buoyancy flux """
function wb(snapshots, i)
    α = parameters.α
    g = parameters.g

    w = snapshots[:w][i]
    T = snapshots[:T][i]

    return α * g * T * w
end

""" horizontal kinetic energy """
function KE(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]

    return 0.5 * (u^2 + v^2)
end

""" mixed layer depth """
function MLD(snapshots, i; threshold = 0.1)
    T = snapshots[:T][i]
    h = MixedLayerDepthOperand(T, abs(threshold))

    return h
end

""" vertical vorticity """
function ζ(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]

    return ∂x(v) - ∂y(u)
end

""" horizontal divergence """
function δ(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]

    return KernelFunctionOperation{Center, Center, Center}(div_xyᶜᶜᶜ, grid, u, v)
end

""" potential vorticity """
function PV(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]
    w = snapshots[:w][i]
    T = snapshots[:T][i]

    f = parameters.f
    α = parameters.α
    g = parameters.g

    ωz = ζ(snapshots, i) + f
    ωx = ∂z(v) - ∂y(w)
    ωy = ∂x(w) - ∂z(u)

    bx = α * g * ∂x(T)
    by = α * g * ∂y(T)
    bz = α * g * ∂z(T)

    return ωx * bx + ωy * by + ωz * bz
end
