using Documenter, FFTW

makedocs(
    modules = [FFTW],
    clean = false,
    sitename = "FFTW.jl",
    pages = Any[
        "Home" => "index.md",
        "API" => "fft.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaMath/FFTW.jl.git",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
    target = "build",
    deps = nothing,
    make = nothing,
)
