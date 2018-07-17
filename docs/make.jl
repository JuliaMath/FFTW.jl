using Documenter, FFTW

makedocs(
    modules = [FFTW],
    clean = false,
    format = :html,
    sitename = "FFTW.jl",
    pages = Any[
        "Home" => "index.md",
        "API" => "fft.md",
    ],
)

deploydocs(
    julia = "nightly",
    repo = "github.com/JuliaMath/FFTW.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
