using Documenter, DocumenterInterLinks, FFTW

links = InterLinks(
    "Julia" => "https://docs.julialang.org/en/v1/",
    "AbstractFFTs" => "https://juliamath.github.io/AbstractFFTs.jl/dev/",
)

makedocs(
    modules = [FFTW],
    clean = false,
    repo = Remotes.GitHub("JuliaMath", "FFTW.jl"),
    sitename = "FFTW.jl",
    pages = Any[
        "Home" => "index.md",
        "API" => "fft.md",
        "Examples" => "examples.md",
    ],
    plugins=[links],
)

deploydocs(
    repo = "github.com/JuliaMath/FFTW.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
