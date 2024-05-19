
@inline transformR(r, p) = 2(p.R - r) / p.R * π * p.Lf - π/2 * (p.Lf - 1)

@inline function uᵢ(x, y, z)
    Lf = parameters.Lf
    R  = 25e3

    return eddy_tangential_velocity(x, y, z, R, Lf, sin)
end

@inline function vᵢ(x, y, z)
    Lf = parameters.Lf
    R  = 25e3

    return eddy_tangential_velocity(x, y, z, R, Lf, cos)
end

@inline function eddy_tangential_velocity(x, y, z, R, Lf, trig)
    # divide into 4 regions
    if x < 50e3 && y < 50e3 # Region 1: warm eddy!
        x′ = x - R
        y′ = y - R
        
        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Lf))
        uθ = warm_eddy_velocity(ξ, z, R, Lf)
        θ  = atan(y′, x′)
        return trig(θ) * uθ
    
    elseif x < 50e3 && y >= 50e3 # Region 2: cold eddy!
        x′ = x - R
        y′ = y - 3R

        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Lf))
        uθ = cold_eddy_velocity(ξ, z, R, Lf)
        θ  = atan(y′, x′)
        return trig(θ) * uθ

    elseif x >= 50e3 && y < 50e3 # Region 3: cold eddy!
        x′ = x - 3R
        y′ = y - R

        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Lf))
        uθ = cold_eddy_velocity(ξ, z, R, Lf)
        θ  = atan(y′, x′)
        return trig(θ) * uθ

    else # Region 4: warm eddy!
        x′ = x - 3R
        y′ = y - 3R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Lf))
        uθ = warm_eddy_velocity(ξ, z, R, Lf)
        θ  = atan(y′, x′)
        return trig(θ) * uθ
    end
end

@inline function Tᵢ(x, y, z)

    Lf = parameters.Lf
    R = 25e3

    # divide into 4 regions
    if x < 50e3 && y < 50e3 # Region 1: warm eddy!
        x′ = x - R
        y′ = y - R
        
        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Lf))
        return warm_eddy(ξ, z)
    
    elseif x < 50e3 && y >= 50e3 # Region 2: cold eddy!
        x′ = x - R
        y′ = y - 3R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Lf))
        return cold_eddy(ξ, z)

    elseif x >= 50e3 && y < 50e3 # Region 3: cold eddy!
        x′ = x - 3R
        y′ = y - R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Lf))
        return cold_eddy(ξ, z)

    else # Region 4: warm eddy!
        x′ = x - 3R
        y′ = y - 3R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Lf))
        return warm_eddy(ξ, z)
    end
end

@inline function T̅⁺(ξ) 
    ΔT = parameters.ΔT 
    T₀ = parameters.T₀
    
    T₃ = 1 - (π - ξ - sin(π - ξ) * cos(π - ξ)) / π
    Tₘ = Int(ξ > 3.1415926535897) + Int(0 < ξ < 3.1415926535897) * T₃

    return ΔT * Tₘ + T₀
end

@inline function T̅⁻(ξ) 
    ΔT = parameters.ΔT 
    T₀ = parameters.T₀
    
    T₃ = 1 - (π - ξ - sin(π - ξ) * cos(π - ξ)) / π
    Tₘ = Int(ξ > 3.1415926535897) + Int(0 < ξ < 3.1415926535897) * T₃

    return T₀ - ΔT * Tₘ
end

# Mixed layer depth profile
@inline function h̅⁺(ξ)
    Δh = parameters.Δm
    h₀ = parameters.m₀

    hᵪ = 1 - (π - ξ - sin(π - ξ) * cos(π - ξ)) / π
    hₘ = Int(ξ > 3.1415926535897) + Int(0 < ξ < 3.1415926535897) * hᵪ

    return Δh * hₘ + h₀
end

@inline function h̅⁻(ξ)
    Δh = parameters.Δm
    h₀ = parameters.m₀

    hᵪ = 1 - (π - ξ - sin(π - ξ) * cos(π - ξ)) / π
    hₘ = Int(ξ > 3.1415926535897) + Int(0 < ξ < 3.1415926535897) * hᵪ

    return h₀ - Δh * hₘ
end

@inline function warm_eddy(ξ, z)

    Lz = parameters.Lz
    T₀ = parameters.T₀
    ΔT = parameters.ΔT

    h = h̅⁺(ξ)
    if z > - h
        return T̅⁺(ξ)
    else
        return (T̅⁺(ξ) - T₀ + 1.2ΔT) / (Lz - h)^2 * (Lz + z)^2 + T₀ - 1.2ΔT
    end
end

@inline function cold_eddy(ξ, z)

    Lz = parameters.Lz
    T₀ = parameters.T₀
    ΔT = parameters.ΔT

    h = h̅⁻(ξ)
    if z > - h
        return T̅⁻(ξ)
    else
        return (T̅⁻(ξ) - T₀ + 1.2ΔT) / (Lz - h)^2 * (Lz + z)^2 + T₀ - 1.2ΔT
    end
end

@inline function warm_eddy_velocity(ξ, z, R, Lf)

    Lz = parameters.Lz
    ΔT = parameters.ΔT
    f  = parameters.f
    α  = parameters.α
    g  = parameters.g

    ∂b∂ξ = - g * α * ΔT * (sin(ξ)^2 - cos(ξ)^2 + 1) / π
    ∂b∂ξ = Int(0 < ξ < 3.1415926535897) * ∂b∂ξ
    ∂ξ∂r = - 2π / R * Lf

    h = h̅⁺(ξ)
    if z > - h
        return ∂ξ∂r * ∂b∂ξ / f * (Lz - h) / 3
    else
        return ∂ξ∂r * ∂b∂ξ / f * (Lz + z)^3 / (Lz - h)^2 / 3
    end
end

@inline function cold_eddy_velocity(ξ, z, R, Lf)

    Lz = parameters.Lz
    ΔT = parameters.ΔT
    f  = parameters.f
    α  = parameters.α
    g  = parameters.g

    ∂b∂ξ = g * α * ΔT * (sin(ξ)^2 - cos(ξ)^2 + 1) / π
    ∂b∂ξ = Int(0 < ξ < 3.1415926535897) * ∂b∂ξ
    ∂ξ∂r = - 2π / R * Lf

    h = h̅⁻(ξ)
    if z > - h
        return ∂ξ∂r * ∂b∂ξ / f * (Lz - h) / 3
    else
        return ∂ξ∂r * ∂b∂ξ / f * (Lz + z)^3 / (Lz - h)^2 / 3
    end
end
