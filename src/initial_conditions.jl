function 𝒴(y) 
    Ly = problem_constants.Ly
    Lf = problem_constants.Lf
    if y ≤ Ly / 2
        return 0.5 * (1 - tanh(y / Lf) + tanh((y - Ly / 2) / Lf))
    else
        return 0.5 * (tanh((y - Ly / 2) / Lf) - tanh((y - Ly) / Lf) - 1)
    end
end

function Φ(x, y, z)
    α  = problem_constants.α
    Lx = problem_constants.Lx
    Ly = problem_constants.Ly
    Lf = problem_constants.Lf

    return α * Lx * Ly / (8π^2) * cos(2π * x / Lx) * sin(4π * (y - Ly / 2) / Ly)
end
