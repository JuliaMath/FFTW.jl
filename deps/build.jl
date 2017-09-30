# If BLAS was compiled with MKL, we want to use MKL for FFTs as well. Thus
# we have to do this little dance to get around having to use BinDeps for
# a library that's already linked to Julia.
if Base.BLAS.vendor() === :mkl
    mklpath = Libdl.dlpath("libmkl_rt")
    depsfile = joinpath(@__DIR__, "deps.jl")
    isfile(depsfile) && rm(depsfile, force=true)
    open(depsfile, "w") do f
        println(f, """
            # This is an auto-generated file, do not edit
            if Libdl.dlopen_e("$mklpath") == C_NULL
                error("Unable to load MKL from '$mklpath'.\\n",
                      "Please rerun Pkg.build(\\"FFTW\\") and restart Julia.")
            end
            const libfftw = "$mklpath"
            const libfftwf = "$mklpath"
        """)
    end
else
    include("build_fftw.jl")
end
