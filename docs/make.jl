using Documenter, FFTW

makedocs(
    modules = [FFTW],
    clean = false,
    format = :html,
    sitename = "FFTW.jl",
    pages = Any[
        "Home" => "index.md",
        "Fourier Transforms" => "fft.md",
        "Signal Processing" => "dsp.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaMath/FFTW.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
