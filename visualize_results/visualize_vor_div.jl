
using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!

using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, ζ, δ
set_theme!(Theme(fontsize = 12))

function visualize_vorticity_divergence(cooling, wind, dTf)
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
    snapshots = load_snapshots(filename; metadata)

    # Let's pick the last snapshot!
    times = snapshots[:T].times
    snapshot_number = length(times)÷2
    nday = @sprintf("%.0f", (times[snapshot_number])/60^2/24)
    println("Plotting snapshot $snapshot_number on day $(nday)...")
    t0 = now()

    # Plot the vorticity and divergence
    f = parameters.f
    # Coordinate arrays
    x, y, z = nodes(snapshots[:T][end])
    
    #xζ, yζ, zζ = nodes(ω)
    #xδ, yδ, zδ = nodes(d)
    ks = 113 # index of the surface layer
    fig = Figure(size = (800, 420))
    axis_kwargs = (xlabel = "x (km)", ylabel = "y (km)",
                  limits = ((0, 100), (0, 100)),aspect = AxisAspect(1))
    axis_kwargs2 = NamedTuple{(:xlabel,:limits,:aspect)}(axis_kwargs)
    ax_ζ = Axis(fig[2, 1][1,1]; title=L"\zeta/f", axis_kwargs...)
    ax_δ = Axis(fig[2, 2][1,1]; title=L"\delta/f", axis_kwargs2...)

    vbnd1,vbnd2 = 7,3
    n = Observable(1)
    ω = @lift interior(compute!(Field(ζ(snapshots, $n))),:,:,ks)/f
    d = @lift interior(compute!(Field(δ(snapshots, $n))),:,:,ks)/f
    #println("Computing vorticity wall time: $((now() - t0).value/1e3) seconds.")
    #println("Computing divergence wall time: $((now() - t0).value/1e3) seconds.")
    hm_ζ = heatmap!(ax_ζ, 1e-3x, 1e-3y, ω; rasterize = true, 
                    colormap = :balance, colorrange = (-vbnd1, vbnd1))
    hm_δ = heatmap!(ax_δ, 1e-3x, 1e-3y, d; rasterize = true, 
                    colormap = :balance, colorrange = (-vbnd2, vbnd2))
    hideydecorations!(ax_δ, ticks = false)
    Colorbar(fig[2, 1][1, 2], hm_ζ)
    Colorbar(fig[2, 2][1, 2], hm_δ)
    
    title = @lift "t = " * string(round(times[$n]/3600/24, digits=2)) * " days"
    Label(fig[1, 1:2], title, fontsize=16, tellwidth=false)
    frames = 1:16:length(times)

    @info "Making a neat animation of vorticity and divergence..."

    record(fig, filesave * "vordiv_" * fileparams * ".mp4", frames, framerate=4) do i
        n[] = i
    end

   # save(filesave * "vordiv_" * fileparams * "_d$(nday).pdf", fig)
    return
end

cooling, wind, dTf = 50, 0.2, 2
visualize_vorticity_divergence(cooling, wind, dTf)