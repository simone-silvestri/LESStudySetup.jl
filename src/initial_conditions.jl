""" horizontal fronts """
function 𝒴(y) 
    Ly = parameters.Ly
    Lf = parameters.Lf
    if y ≤ Ly / 2
        return 0.5 * (1 - tanh(y / Lf) + tanh((y - Ly / 2) / Lf))
    else
        return 0.5 * (tanh((y - Ly / 2) / Lf) - tanh((y - Ly) / Lf) - 1)
    end
end

""" initial barotropic streamfunction """
function Φ(x, y, z)
    α  = parameters.α
    Lx = parameters.Lx
    Ly = parameters.Ly
    Lf = parameters.Lf

    return α * Lx * Ly / (8π^2) * cos(2π * x / Lx) * sin(4π * (y - Ly / 2) / Ly)
end

""" initial buoyancy field """
function bᵢ(x, y, z)
    Δρ = parameters.Δρ
    ρ₀ = parameters.ρ₀ 
    N² = parameters.N² 

    ρ′ = Δρ * 𝒴(y) - ρ₀ / 9.80655 * N² * z

    return - 9.8655 * ρ′ / ρ₀
end

""" initial temperature field Tᵢ = bᵢ / (α ⋅ g) """
Tᵢ(x, y, z) = bᵢ(x, y, z) / parameters.α / 9.80655 + 19

""" initial zonal velocity uᵢ = - ∂yΦ """
function uᵢ(x, y, z)
    θ  = parameters.θ
    Lx = parameters.Lx
    Ly = parameters.Ly

    return - θ * Lx * Ly / (8π^2) * cos(2π * x / Lx) * 4π / Ly * cos(4π * (y - Ly / 2) / Ly)
end

""" initial meridional veloctity vᵢ = ∂xΦ """
function vᵢ(x, y, z)
    θ  = parameters.θ
    Lx = parameters.Lx
    Ly = parameters.Ly

    return - θ * Lx * Ly / (8π^2) * 2π / Lx * sin(2π * x / Lx) * sin(4π * (y - Ly / 2) / Ly)
end

