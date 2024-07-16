using FFTW
using Oceananigans.Grids: φnode

struct Spectrum{S, F}
    spec :: S
    freq :: F
end

import Base

Base.:(+)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .+ t.spec, s.freq)
Base.:(*)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .* t.spec, s.freq)
Base.:(/)(s::Spectrum, t::Int)      = Spectrum(s.spec ./ t, s.freq)

Base.real(s::Spectrum) = Spectrum(real.(s.spec), s.freq)
Base.abs(s::Spectrum)  = Spectrum( abs.(s.spec), s.freq)

function power_cospectrum_1d(var1, var2, x)

    Nx  = length(x)
    Nfx = Int64(Nx)
    
    spectra = zeros(ComplexF64, Int(Nfx/2))
    
    dx = x[2] - x[1]

    freqs = fftfreq(Nfx, 1.0 / dx) # 0, +ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfx/2)] .* 2.0 .* π
    
    fourier1   = fft(var1) / Nfx
    fourier2   = fft(var2) / Nfx
    spectra[1] += fourier1[1] .* conj(fourier2[1]) .+ fourier2[1] .* conj(fourier1[1])

    for m in 2:Int(Nfx/2)
        spectra[m] += fourier1[m] .* conj(fourier2[m]) .+ fourier2[m] .* conj(fourier1[m])
    end
    return Spectrum(spectra, freqs)
end

