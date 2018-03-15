using Libdl

# If BLAS was compiled with MKL and the user wants MKL-based FFTs, we'll oblige.
# In that case, we have to do this little dance to get around having to use BinDeps
# for a library that's already linked to Julia.
settings = joinpath(@__DIR__, "..", ".build_settings")
if haskey(ENV, "JULIA_FFTW_PROVIDER")
    provider = ENV["JULIA_FFTW_PROVIDER"]
    open(f -> println(f, provider), settings, "w")
elseif isfile(settings)
    provider = readchomp(settings)
else
    provider = "FFTW"
    open(f -> println(f, provider), settings, "w")
end
if provider == "MKL" && Base.BLAS.vendor() === :mkl
    mklpath = Libdl.dlpath("libmkl_rt")
    depsfile = joinpath(@__DIR__, "deps.jl")
    isfile(depsfile) && rm(depsfile, force=true)
    open(depsfile, "w") do f
        println(f, """
            # This is an auto-generated file, do not edit
            using Libdl
            if Libdl.dlopen_e("$mklpath") == C_NULL
                error("Unable to load MKL from '$mklpath'.\\n",
                      "Please rerun Pkg.build(\\"FFTW\\") and restart Julia.")
            end
            const libfftw = "$mklpath"
            const libfftwf = "$mklpath"
        """)
    end
elseif provider == "MKL"
    error("MKL build requested for FFTW but Julia was not built with MKL.\n",
          "To fix this, set ENV[\"JULIA_FFTW_PROVIDER\"] = \"FFTW\" and \n",
          "rerun Pkg.build(\"FFTW\").")
else
    include("build_fftw.jl")
end
