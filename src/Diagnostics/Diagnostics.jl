module Diagnostics

using Oceananigans
using Oceananigans

using Oceananigans
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Architectures: device, architecture
using Oceananigans.Utils: launch!
using Oceananigans.Grids: Center, Face, inactive_node, znode
using Oceananigans.Operators: Δzᶜᶜᶠ, ζ₃ᶠᶠᶜ

using KernelAbstractions: @index, @kernel
using KernelAbstractions.Extras.LoopInfo: @unroll

import Oceananigans.Fields: compute!

using Oceananigans.Fields: OneField, condition_operand
using Oceananigans.AbstractOperations: materialize_condition!
using Oceananigans.Utils

function load_snapshots(filename; architecture = CPU())

    snapshots = Dict()

    u = FieldTimeSeries(filename, "u"; architecture, backend = OnDisk())
    v = FieldTimeSeries(filename, "v"; architecture, backend = OnDisk())
    w = FieldTimeSeries(filename, "w"; architecture, backend = OnDisk())
    T = FieldTimeSeries(filename, "T"; architecture, backend = OnDisk())

    snapshots[:u] = u
    snapshots[:v] = v
    snapshots[:w] = w
    snapshots[:T] = T
    
    return snapshots
end

function VerticalVorticity(u, v)
    ζ = Field(∂x(v) - ∂y(u))
    compute!(ζ)
    return ζ
end

include("mixed_layer.jl")

end