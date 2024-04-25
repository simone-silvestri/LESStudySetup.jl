####
#### Double temperature front
####

@inline transformX(x, p) = x ≤ p.Lx / 2 ? 
                        2x / p.Lx * π * (1 + p.Lf) - π/2 * p.Lf :
                        2(p.Lx - x) / p.Lx * π * (1 + p.Lf) - π/2 * p.Lf

""" background temperature """
@inline function T̅(x, y, z, p) 
    ΔT = p.ΔT 
    T₀ = p.T₀
    
    ξ = transformX(x, p)
    T = ifelse(ξ < 0, 0, ifelse(ξ > π, 1, 1 - (π - ξ - sin(π - ξ) * cos(π - ξ)) / π))
    return ΔT * T + T₀
end

norm_x(x, p) = 2 * (x - p.Lx / 2) / p.Lx # from -1 to 1
norm_y(y, p) = 2 * (y - p.Ly / 2) / p.Ly # from -1 to 1

@inline η(x, y, p) = exp(-(norm_x(x, p)^2 + norm_y(y, p)^2) ./ p.σ²) 

""" background and initial zonal velocity """
@inline function U̅(x, y, z, p)
    f = p.f
    g = p.g
    return - g * 2 * norm_y(y, p) * η(x, y, p) / f / p.Ly * 2
end

""" background and initial meridional velocity """
@inline function V̅(x, y, z, p)
    f  = p.f
    ΔT = p.ΔT
    g  = p.g
    
    ξ = transformX(x, p)
    ∂b∂ξ = ifelse(ξ < 0, 0, ifelse(ξ > π, 0, (sin(ξ)^2 - cos(ξ)^2 + 1) / π))
    ∂ξ∂x = 2π / p.Lx

    return g * 2 * norm_x(x, p) * η(x, y, p) / f / p.Lx * 2 + ΔT * ∂b∂ξ * ∂ξ∂x * (p.Lz + z)
end

""" initial temperature field """
@inline function Tᵢ(x, y, z, p)

    N² = p.N²
    α  = p.α
    g  = p.g
    H  = p.H
    ΔH = p.ΔH

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
