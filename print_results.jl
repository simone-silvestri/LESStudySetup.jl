using CairoMakie, JLD2

file = jldopen("hydrostatic_snapshots_expweak_cooling_weak_wind_mixed.jld2") 

iterations = parse.(Int, keys(file["timeseries/t"]))
iter = Observable(iterations[1]);
ui = @lift(file["timeseries/u/" * string($iter)][:, :, 120])
vi = @lift(file["timeseries/v/" * string($iter)][:, :, 120])
Ti = @lift(file["timeseries/T/" * string($iter)][:, :, 120])
wi = @lift(file["timeseries/w/" * string($iter)][:, :, 70])

ζi = @lift begin
       um = $ui[:, 2:end] .- $ui[:, 1:end-1]
       vm = $vi[2:end, :] .- $vi[1:end-1, :]
       um = um[1:end-1, :]
       vm = vm[:, 1:end-1]
       (vm .- um) / 200
end

fig = Figure(size = (1400, 300))
ax = Axis(fig[1, 1], title = "Vertical velocity");
heatmap!(ax, wi, colorrange = (-1e-2, 1e-2))
ax = Axis(fig[1, 2], title = "Vertical vorticity");
heatmap!(ax, ζi, colorrange = (-5e-4, 5e-4), colormap = :bwr)
ax = Axis(fig[1, 3], title = "Temperature");
heatmap!(ax, Ti, colorrange = (4.8, 6.5), colormap = :magma)

CairoMakie.record(fig, "weak_forcing.mp4", iterations; framerate = 8) do i
       @info "doing iter $i"
       iter[] = i
end
