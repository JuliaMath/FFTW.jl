__precompile__()

module FFTW

const depsfile = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
if isfile(depsfile)
    include(depsfile)
else
    error("FFTW is not properly installed. Please run Pkg.build(\"FFTW\") ",
          "and restart Julia.")
end

# MKL provides its own FFTW
fftw_vendor() = Base.BLAS.vendor() === :mkl ? :mkl : :fftw

if fftw_vendor() === :mkl
    const libfftw_name = "libmkl_rt"
    const libfftwf_name = "libmkl_rt"
elseif is_windows()
    const libfftw_name = "libfftw3"
    const libfftwf_name = "libfftw3f"
else
    const libfftw_name = "libfftw3_threads"
    const libfftwf_name = "libfftw3f_threads"
end

include("dft.jl")
include("fft.jl")
include("dct.jl")
include("dsp.jl") # TODO: Move these functions to DSP.jl

end # module
