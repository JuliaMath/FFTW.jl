using LinearAlgebra

import Libdl
const depsfile = joinpath(@__DIR__, "deps.jl")

settings = joinpath(first(DEPOT_PATH), "prefs", "FFTW")
mkpath(dirname(settings))
if haskey(ENV, "JULIA_FFTW_PROVIDER")
    provider = ENV["JULIA_FFTW_PROVIDER"]
    open(f -> println(f, provider), settings, "w")
elseif isfile(settings)
    provider = readchomp(settings)
else
    provider = "FFTW"
    open(f -> println(f, provider), settings, "w")
end

if provider == "MKL"
    const mkllib = Sys.iswindows() ? "mkl_rt" : "libmkl_rt"
    # If BLAS was compiled with MKL and the user wants MKL-based FFTs, we'll oblige.
    open(depsfile, "w") do io
        println(io, """
            using IntelOpenMP_jll, MKL_jll
            check_deps() = nothing
            const libfftw3 = MKL_jll.libmkl_rt_path
            const libfftw3f = libfftw3
        """)
    end
elseif provider == "FFTW"
    open(depsfile, "w") do io
        println(io, """
            using FFTW_jll
            check_deps() = nothing
        """)
    end
else
    error("Unrecognized JULIA_FFTW_PROVIDER \"$provider\".\n",
          "To fix this, set ENV[\"JULIA_FFTW_PROVIDER\"] to \"FFTW\" or \"MKL\"\n",
          "and rerun Pkg.build(\"FFTW\").")
end
