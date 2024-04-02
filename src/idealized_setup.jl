
function idealized_setup(arch; hydrostatic_approximation = true)
    
    # Retrieving the problem constants
    Δbʸ = problem_constants.Δbʸ
    ρ₀  = problem_constants.ρ₀ 
    N²  = problem_constants.N² 
    Δh  = problem_constants.Δh 
    Lx  = problem_constants.Lx 
    Ly  = problem_constants.Ly 
    Lz  = problem_constants.Lz 
    Lf  = problem_constants.Lf 
     α  = problem_constants.α  
     f  = problem_constants.f  
    τw  = problem_constants.τw 

    Nx = ceil(Int, Lx / Δh)
    Ny = ceil(Int, Lx / Δh)
    Nz = ceil(Int, Lx / Δh)


end