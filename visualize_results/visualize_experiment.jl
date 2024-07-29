using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, spatial_filtering, power_cospectrum_1d, isotropic_powerspectrum, ζ, δ
set_theme!(Theme(fontsize = 12))

function visualize(cooling, wind, dTf)
    # Examples! (fill in the correct filename and metadata filename)
    # cooling, wind, dTf = 25, 0.02, -1
    cooling = @sprintf("%03d", cooling)
    wind = replace("$(wind)","." => "" )
    if dTf < 0
        fileparams = "four_vortices_cooling_$(cooling)_wind_$(wind)"
    else
        if length(wind) < 2
            wind = "0" * wind
        end
        fileparams = "free_surface_short_test_$(cooling)_wind_$(wind)_dTf_$(dTf)"
    end
    filehead = "/orcd/data/abodner/001/simone/LESStudySetup.jl/experiments/"
    filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
    metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
    filesave = "/orcd/data/abodner/001/shirui/LESStudySetup/results/"

    # load all the data!!
    println("Loading data from $filename...")
    snapshots = load_snapshots(filename; metadata)

    # Let's pick the last snapshot!
    times = snapshots[:T].times
    snapshot_number = length(times) ÷ 2
    nday = @sprintf("%2.0f", (times[snapshot_number])/60^2/24)
    println("Plotting snapshot $snapshot_number on day $(nday)...")
    t0 = now()

    T = snapshots[:T][snapshot_number]
    u = snapshots[:u][snapshot_number]
    v = snapshots[:v][snapshot_number]
    w = snapshots[:w][snapshot_number]
    #s = compute!(Field(sqrt(u^2 + v^2)))
    println("Loading fields wall time: $((now() - t0).value/1e3) seconds.")

    grid = T.grid

    # Coordinate arrays
    xu, yu, zu = nodes(u)
    xv, yv, zv = nodes(v)
    #xs, ys, zs = nodes(s)
    xw, yw, zw = nodes(w)
    xT, yT, zT = nodes(T)

    # Plot the fields
    kT, ks, kw = 113, 124, 113
    jslices = [250,500,750]
    Tmin, Tmax = minimum(interior(T)), maximum(interior(T))
    vbnd, wbnd = maximum(abs, interior(v)), 0.005

    fig = Figure(size = (800, 250))
    axis_kwargs = (xlabel = "x (km)", ylabel = "y (km)",
                limits = ((0, 100), (0, 100)))
    axis_kwargs2 = NamedTuple{(:xlabel,:limits)}(axis_kwargs)
    ax_T = Axis(fig[1, 1][1,1]; title="Temperature (ᵒC), z=$(zT[kT])m", axis_kwargs...)
    #ax_s = Axis(fig[1, 2][1,1]; title="Horizontal speed (m s⁻¹)", axis_kwargs2...)
    ax_v = Axis(fig[1, 2][1,1]; title="Meridional velocity (m s⁻¹), z=$(zv[ks])m", axis_kwargs2...)
    ax_w = Axis(fig[1, 3][1,1]; title="Vertical velocity (m s⁻¹), z=$(zw[kw])m", axis_kwargs2...)    
    hm_T = heatmap!(ax_T, 1e-3xT, 1e-3yT, interior(T,:,:,kT); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
    Colorbar(fig[1, 1][1, 2], hm_T)
    #= hm_s = heatmap!(ax_s, 1e-3xs, 1e-3ys, interior(s,:,:,ks); colormap = :speed, colorrange = (0, 0.85))
    hideydecorations!(ax_s, ticks = false)
    Colorbar(fig[1, 2][1, 2], hm_s) =#
    hm_v = heatmap!(ax_v, 1e-3xv, 1e-3yv, interior(v,:,:,ks); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
    hideydecorations!(ax_v, ticks = false)
    Colorbar(fig[1, 2][1, 2], hm_v)
    hm_w = heatmap!(ax_w, 1e-3xw, 1e-3yw, interior(w,:,:,kw); rasterize = true, colormap = :balance, colorrange = (-wbnd, wbnd))#, highclip = :red, lowclip = :blue)
    hideydecorations!(ax_w, ticks = false)
    Colorbar(fig[1, 3][1, 2], hm_w)
    hlines!(ax_T, 1e-3*yT[jslices]; color = :black, linewidth = 0.5)
    hlines!(ax_v, 1e-3*yv[jslices]; color = :black, linewidth = 0.5)
    hlines!(ax_w, 1e-3*yw[jslices]; color = :black, linewidth = 0.5)
    save(filesave * "fields_" * fileparams * "_d$(nday).pdf", fig)
    println("Finished plotting fields, wall time: $((now() - t0).value/1e3) seconds.")
    zmin = -250
    axis_kwargs0 = (xlabel = "x (km)",
                    ylabel = "z (m)",
                    limits = ((0, 100), (zmin, 0)))
    axis_kwargs1 = NamedTuple{(:xlabel,:limits)}(axis_kwargs0)
    axis_kwargs2 = NamedTuple{(:ylabel,:limits)}(axis_kwargs0)
    #axis_kwargs3 = NamedTuple{(:limits)}(axis_kwargs0)
    α = parameters.α
    ρ₀, Δρ = 1020, 0.03
    fig = Figure(size = (800, 400))
    for j in [3,2,1]
        jT, jv, jw = 250*j, 250*j, 250*j
        if j == 3
            local ax_T = Axis(fig[j, 1][1,1]; axis_kwargs0...)
            local ax_v = Axis(fig[j, 2][1,1]; axis_kwargs1...)
            local ax_w = Axis(fig[j, 3][1,1]; axis_kwargs1...)
        elseif j == 2
            local ax_T = Axis(fig[j, 1][1,1]; axis_kwargs2...)
            local ax_v = Axis(fig[j, 2][1,1]; limits = ((0, 100), (zmin, 0)))
            local ax_w = Axis(fig[j, 3][1,1]; limits = ((0, 100), (zmin, 0)))
        else
            local ax_T = Axis(fig[j, 1][1,1]; title="Temperature (ᵒC)", axis_kwargs2...)
            local ax_v = Axis(fig[j, 2][1,1]; title="Meridional velocity (m s⁻¹)", limits = ((0, 100), (zmin, 0)))
            local ax_w = Axis(fig[j, 3][1,1]; title="Vertical velocity (m s⁻¹)", limits = ((0, 100), (zmin, 0)))
        end
        
        hm_T = heatmap!(ax_T, 1e-3xT, zT, interior(T,:,jT,:); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
        hm_v = heatmap!(ax_v, 1e-3xv, zv, interior(v,:,jv,:); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
        hm_w = heatmap!(ax_w, 1e-3xw, zw, interior(w,:,jw,:); rasterize = true, colormap = :balance, colorrange = (-wbnd, wbnd))#, highclip = :red, lowclip = :blue)
        if j < 3
            hidexdecorations!(ax_T, ticks = false)
            hidexdecorations!(ax_v, ticks = false)
            hidexdecorations!(ax_w, ticks = false)
        end
        hideydecorations!(ax_v, ticks = false)
        hideydecorations!(ax_w, ticks = false)
        Colorbar(fig[j, 1][1, 2], hm_T)
        Colorbar(fig[j, 2][1, 2], hm_v)
        Colorbar(fig[j, 3][1, 2], hm_w)
        hlines!(ax_T, [zT[kT]]; color = :black, linewidth = 0.5)
        hlines!(ax_v, [zv[ks]]; color = :black, linewidth = 0.5)
        hlines!(ax_w, [zw[kw]]; color = :black, linewidth = 0.5)

        ρj = -ρ₀*α*interior(T,:,jT,:)
        ρj10 = (ρj[:,120].+ρj[:,121])/2
        h = zeros(grid.Nx)
        for i in 1:grid.Nx
            ki = findfirst(ρj[i,1:119] .<= ρj10[i]+Δρ)
            si = (zT[ki]-zT[ki-1])/(ρj[i,ki]-ρj[i,ki-1])
            h[i] = -(zT[ki] + si*(ρj10[i]+Δρ-ρj[i,ki]))
        end
        println("mean mixed layer depth along y = $(1e-3yT[jT])km is $(mean(h))m")
        lines!(ax_T, -h; color = :white, linewidth = 0.8)
    end
    save(filesave * "fieldslices_" * fileparams * "_d$(nday).pdf", fig)
    println("Finished plotting fieldslices, wall time: $((now() - t0).value/1e3) seconds.")

    return
end

coolings = [50,50,0]
winds = [0.1,0.2,0]
dTfs = [2,2,2]
for i in 1:length(coolings)
    visualize(coolings[i], winds[i], dTfs[i])
end

#= # Filtering the fields to find mean values
T̅ = spatial_filtering(T)
U = spatial_filtering(u)
V = spatial_filtering(v)
W = spatial_filtering(w)

# Calculate the fluctuations
T′ = compute!(Field(T - T̅))
u′ = compute!(Field(u - U))
v′ = compute!(Field(v - V))
w′ = compute!(Field(w - W))

# Compute eddy fluxes 
u′T′ = compute!(Field(u′ * T′)) 
v′T′ = compute!(Field(v′ * T′)) 
w′T′ = compute!(Field(w′ * T′))  =#

