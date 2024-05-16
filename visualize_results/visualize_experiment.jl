using LESStudySetup
using CairoMakie
using SixelTerm

using Statistics: mean
using Oceananigans: compute!
using LESStudySetup

using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots

# # Retrieve the grid from the snapshots
# grid = T.grid 
# Nx, Ny, Nz = size(grid)

# # How many snapshots are there?
# number_of_snapshots = length(u.times)

# iter = Observable(1)

# # Produce a nice video of surface temperature evolution
# Ti = @lift(interior(T[$iter], :, :, grid.Nz)) # Surface field!
# fig = Figure()
# ax = Axis(fig[1, 1])
# heatmap!(ax, Ti)
# CairoMakie.record(fig, "center-eddy-T.mp4", 1:number_of_snapshots, framerate = 5) do i
#     @info "iter $i"
#     iter[] = i
# end

# # Calulate the mixed layer depth at the last timestep
# T_end = T[number_of_snapshots];
# mixed_layer = MixedLayerDepth(grid, (; T = T_end); ΔT = 0.00002)
# compute!(mixed_layer)

# # Check the meridionally average temperature profile
# Tmean = @lift(mean(interior(T[$iter]), dims = 2)[:, 1, :])
# fig = Figure()
# ax = Axis(fig[1, 1])
# heatmap!(ax, Tmean)
# CairoMakie.record(fig, "center-eddy-T.mp4", 1:number_of_snapshots, framerate = 5) do i
#     @info "iter $i"
#     iter[] = i
# end

# # Check the vorticity
# ζ = @lift begin
#     z = Diagnostics.VerticalVorticity(u[$iter], v[$iter])
#     interior(z, :, :, 125)
# end

# fig = Figure()
# ax = Axis(fig[1, 1])
# heatmap!(ax, ζ, colormap = :bwr, clims = (-1e-5, 1e-5)
# CairoMakie.record(fig, "vorticity.mp4", 1:100, framerate = 5) do i
#     @info "iter $i"
#     iter[] = i
# end


