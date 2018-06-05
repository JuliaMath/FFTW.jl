import Libdl
const depsfile = joinpath(@__DIR__, "deps.jl")

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
    mklpath = escape_string(Libdl.dlpath("libmkl_rt"))
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
elseif provider == "MKL"
    error("MKL build requested for FFTW but Julia was not built with MKL.\n",
          "To fix this, set ENV[\"JULIA_FFTW_PROVIDER\"] = \"FFTW\" and \n",
          "rerun Pkg.build(\"FFTW\").")
elseif provider != "FFTW"
    error("Unrecognized JULIA_FFTW_PROVIDER \"$provider\".\n",
          "To fix this, set ENV[\"JULIA_FFTW_PROVIDER\"] to \"FFTW\" or \"MKL\"\n",
          "and rerun Pkg.build(\"FFTW\").")
else
    include("build_fftw.jl")
end

# on some platforms, FFTW may be built without threads.  We need to define
# a variable for this in deps.jl so that we can check for it statically,
# in order to decide whether to ccall the FFTW threads functions.
include(depsfile)
open(depsfile, "a") do f
    has_threads = Libdl.dlsym_e(Libdl.dlopen(libfftw3), "fftw_init_threads") != C_NULL
    println(f, "const has_threads = $has_threads")
end
