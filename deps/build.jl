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
    if BLAS.vendor() === :mkl
        mklpath = Libdl.dlpath(mkllib)
    else
        using Conda
        Conda.add("mkl_fft")
        mklpath = joinpath(Conda.lib_dir(Conda.ROOTENV), mkllib)
    end
    mklpath = escape_string(mklpath)
    isfile(depsfile) && rm(depsfile, force=true)
    open(depsfile, "w") do f
        println(f, """
            # This is an auto-generated file, do not edit
            import Libdl
            const libfftw3 = "$mklpath"
            const libfftw3f = libfftw3
            function check_deps()
                if Libdl.dlopen_e(libfftw3) == C_NULL
                    error("Unable to load MKL from '$mklpath'.\\n",
                          "Please rerun Pkg.build(\\"FFTW\\") and restart Julia.")
                end
            end
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
