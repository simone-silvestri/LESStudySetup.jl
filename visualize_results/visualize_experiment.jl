using LESStudySetup
using CairoMakie
using SixelTerm

using Statistics: mean
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using LESStudySetup

using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, spatial_filtering, power_cospectrum_1d

# Examples! (fill in the correct filename and metadata filename)
filename = "hydrostatic_snapshots_free_surface_short_test_050_wind_01_dTf_2.jld2"
metadata = "experiment_free_surface_short_test_050_wind_01_dTf_2_metadata.jld2"

# load all the data!!
snapshots = load_snapshots(filename; metadata)

# Let's pick the last snapshot!
times = snapshots[:T].times
snapshot_number = length(times)

T = snapshots[:T][snapshot_number]
u = snapshots[:u][snapshot_number]
v = snapshots[:v][snapshot_number]
w = snapshots[:w][snapshot_number]

grid = T.grid

# Filtering the fields to find mean values
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
w′T′ = compute!(Field(w′ * T′)) 

# Compute the horizontal spectrum of KE -> [0.5 ⋅ (ûû⋆ + v̂v̂⋆)] at k = 100
klev = 100

KE = power_cospectrum_1d(interior(u, :, 1, klev), interior(u, :, 1, klev), xnodes(u)) / grid.Ny

for j in 2:grid.Ny
    KE += power_cospectrum_1d(interior(u, :, j, klev), interior(u, :, j, klev), xnodes(u)) / grid.Ny
end


WT = power_cospectrum_1d(interior(w, :, 1, klev), interior(T, :, 1, klev), xnodes(u)) / grid.Ny

for j in 2:grid.Ny
    KE += power_cospectrum_1d(interior(w, :, j, klev), interior(T, :, j, klev), xnodes(u)) / grid.Ny
end

