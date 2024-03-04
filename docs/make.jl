using Documenter, FFTW

makedocs(
    modules = [FFTW],
    clean = false,
    sitename = "FFTW.jl",
    pages = Any[
        "Home" => "index.md",
        "API" => "fft.md",
        "Examples" => "examples.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaMath/FFTW.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
