####
#### Double temperature front
####

@inline transformX(x, y) = x ≤ parameters.Lx / 2 ? 
                          2x / parameters.Lx * π * (1 + parameters.Lf) - π/2 * parameters.Lf :
                          2(parameters.Lx - x) / parameters.Lx * π * (1 + parameters.Lf) - π/2 * parameters.Lf

""" background temperature """
@inline function T̅(x, y, z) 
    ΔT = parameters.ΔT 
    T₀ = parameters.T₀
    ξ = transformX(x, y)
    T = ifelse(ξ < 0, 0, ifelse(ξ > π, 1, 1 - (π - ξ - sin(π - ξ) * cos(π - ξ)) / π))
    return ΔT * T + T₀
end

norm_x(x) = 2 * (x - parameters.Lx / 2) / parameters.Lx # from -1 to 1
norm_y(y) = 2 * (y - parameters.Ly / 2) / parameters.Ly # from -1 to 1

@inline η(x, y) = exp(-(norm_x(x)^2 + norm_y(y)^2) ./ parameters.σ²) 

@inline function U̅(x, y, z)
    f = parameters.f
    g = parameters.g
    return - g * 2 * norm_y(y) * η(x, y) / f / parameters.Ly * 2
end

@inline function V̅(x, y, z)
    f  = parameters.f
    ΔT = parameters.ΔT
    g  = parameters.g
    
    ξ = transformX(x, y)
    ∂b∂ξ = ifelse(ξ < 0, 0, ifelse(ξ > π, 0, (sin(ξ)^2 - cos(ξ)^2 + 1) / π))
    ∂ξ∂x = 2π / parameters.Lx

    return g * 2 * norm_x(x) * η(x, y) / f / parameters.Lx * 2 + ΔT * ∂b∂ξ * ∂ξ∂x * (parameters.Lz + z)
end

""" initial temperature field """
function Tᵢ(x, y, z)

    N² = parameters.N²
    α  = parameters.α
    g  = parameters.g
    H  = parameters.H
    ΔH = parameters.ΔH

    T_surface = T̅(x, y, z)

    ## Noise with 8 m decay scale
    Ξ(z) = rand() * exp(z / 8)

    dTdz_thermocline   = N² * 10 / (α * g)
    dTdz               = N² / (α * g)
        
    if z ≥ - H
        return T_surface + dTdz * z
    elseif - H - ΔH ≤ z ≤ - H
        return T_surface - dTdz * H + dTdz_thermocline * (z + H)
    else
        return T_surface - dTdz_thermocline * ΔH + dTdz * (z + ΔH)
    end
end
