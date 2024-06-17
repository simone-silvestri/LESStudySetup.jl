using FFTW
using Oceananigans.Utils
using Oceananigans.Units
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index

function spectral_averaging(u::Field; xcutoff = 20kilometer, ycutoff = 20kilometer)

    Δx = u.grid.Δxᶜᵃᵃ
    Δy = u.grid.Δyᵃᶜᵃ
    
    nx  = ceil(Int, xcutoff / Δx)
    ny  = ceil(Int, ycutoff / Δy)

    nx2 = ceil(Int, xcutoff / Δx * π / 2) 
    ny2 = ceil(Int, xcutoff / Δx * π / 2) 

    ū = deepcopy(u)
    d = interior(u)

    d̂ = fft(d)

    @show nx2, ny2, nx, ny

    Nx, Ny, _ = size(u)

    arch = architecture(u)
    grid = u.grid

    launch!(arch, grid, :xyz, _corse_grain!, d̂, nx, ny, nx2, ny2, Nx, Ny)

    d = ifft(d̂)
    
    set!(ū, real.(d))
    fill_halo_regions!(ū)

    return ū
end

@kernel function _corse_grain!(d̂, nx, ny, nx2, ny2, Nx, Ny)
    i, j, k = @index(Global, NTuple)
    i′ = ifelse(i > Nx ÷ 2, i - Nx ÷ 2, Nx ÷ 2 - i)
    j′ = ifelse(j > Ny ÷ 2, j - Ny ÷ 2, Ny ÷ 2 - j)

    # remove the outside frame till nx and ny
    outside_frame   = (i′ <= nx) | (j′ <= ny)
    smooth_region_x = (i′ > nx)  & (i′ <= nx2)
    smooth_region_y = (j′ > ny)  & (j′ <= ny2)

    only_smooth_x = smooth_region_x & !smooth_region_y
    only_smooth_y = smooth_region_y & !smooth_region_x
    smooth_both   = smooth_region_x & smooth_region_y

    scaling_x = (i′ - nx) / (nx2 - nx)
    scaling_y = (j′ - ny) / (ny2 - ny)
    scaling_b = scaling_x * scaling_y
    
    scaling = ifelse(outside_frame, 0, 
              ifelse(only_smooth_x, scaling_x,
              ifelse(only_smooth_y, scaling_y, 
              ifelse(smooth_both, scaling_b, 1))))

    @inbounds d̂[i, j, k] = scaling * d̂[i, j, k]
end
