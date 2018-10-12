using LinearAlgebra

import Libdl
const depsfile = joinpath(@__DIR__, "deps.jl")

# If BLAS was compiled with MKL and the user wants MKL-based FFTs, we'll oblige.
# In that case, we have to do this little dance to get around having to use BinDeps
# for a library that's already linked to Julia.
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
    if BLAS.vendor() === :mkl
        mklpath = Libdl.dlpath("libmkl_rt")
    else
        using Conda
        Conda.add("mkl_fft")
        mklpath = joinpath(Conda.lib_dir(Conda.ROOTENV), "libmkl_rt")
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
elseif provider != "FFTW"
    error("Unrecognized JULIA_FFTW_PROVIDER \"$provider\".\n",
          "To fix this, set ENV[\"JULIA_FFTW_PROVIDER\"] to \"FFTW\" or \"MKL\"\n",
          "and rerun Pkg.build(\"FFTW\").")
else
    include("build_fftw.jl")
end
