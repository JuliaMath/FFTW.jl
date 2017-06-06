__precompile__()

module FFTW

# Since nothing is exported from AbstractFFTs as long as the FFT functionality is
# defined (or deprecated) in Base, we need to be very explicit about the things we
# want to import
import AbstractFFTs: Plan, ScaledPlan,
                     fft, ifft, bfft, fft!, ifft!, bfft!,
                     plan_fft, plan_ifft, plan_bfft, plan_fft!, plan_ifft!, plan_bfft!,
                     rfft, irfft, brfft, plan_rfft, plan_irfft, plan_brfft,
                     fftshift, ifftshift,
                     rfft_output_size, brfft_output_size,
                     plan_inv, normalization

if !isdefined(Base, :FFTW)
    export dct, idct, dct!, idct!, plan_dct, plan_idct, plan_dct!, plan_idct!
end
if !isdefined(Base, :DSP)
    export filt, filt!, deconv, conv, conv2, xcorr
end

const depsfile = joinpath(dirname(@__DIR__), "deps", "deps.jl")
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

include("fft.jl")
include("dct.jl")
include("dsp.jl") # TODO: Move these functions to DSP.jl

end # module
