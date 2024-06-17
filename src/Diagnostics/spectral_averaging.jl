using FFTW
using Oceananigans.Utils
using Oceananigans.Units
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index

function spectral_averaging(u::Field; xcutoff = 8kilometer, ycutoff = 8kilometer)

    Δx = u.grid.Δxᶜᵃᵃ
    Δy = u.grid.Δyᵃᶜᵃ
    
    nx = ceil(Int, xcutoff / Δx)
    ny = ceil(Int, ycutoff / Δy)

    ū = deepcopy(u)
    d = interior(u)

    d̂ = fft(d)

    Nx, Ny, Nz = size(u)

    Nfx  = size(u, 1) ÷ 2 - nx
    Nfy  = size(u, 2) ÷ 2 - ny
    arch = architecture(u)
    grid = u.grid

    xparams = KernelParameters((2nx, Ny, Nz), (Nfx, 0, 0)) 
    yparams = KernelParameters((Nx, 2ny, Nz), (0, Nfy, 0)) 

    launch!(arch, grid, xparams, _corse_grain!, d̂)
    launch!(arch, grid, yparams, _corse_grain!, d̂)

    d = ifft(d̂)
    
    set!(ū, real.(d))
    fill_halo_regions!(ū)

    return ū
end

@kernel function _corse_grain!(d̂)
    i, j, k = @index(Global, NTuple)
    @inbounds d̂[i, j, k] = 0
end

