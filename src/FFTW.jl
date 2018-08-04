module FFTW

using Compat, LinearAlgebra, Reexport
import Libdl
@reexport using AbstractFFTs

import AbstractFFTs: Plan, ScaledPlan,
                     fft, ifft, bfft, fft!, ifft!, bfft!,
                     plan_fft, plan_ifft, plan_bfft, plan_fft!, plan_ifft!, plan_bfft!,
                     rfft, irfft, brfft, plan_rfft, plan_irfft, plan_brfft,
                     fftshift, ifftshift,
                     rfft_output_size, brfft_output_size,
                     plan_inv, normalization

export dct, idct, dct!, idct!, plan_dct, plan_idct, plan_dct!, plan_idct!

const depsfile = joinpath(dirname(@__DIR__), "deps", "deps.jl")
if isfile(depsfile)
    include(depsfile)
else
    error("FFTW is not properly installed. Please run Pkg.build(\"FFTW\") ",
          "and restart Julia.")
end

# MKL provides its own FFTW
const fftw_vendor = occursin("libmkl_rt", libfftw3) ? :mkl : :fftw

# If FFTW was built with threads, then they must be initialized before any FFTW planning routine.
#   -- This initializes FFTW's threads support (defaulting to 1 thread).
#      If this isn't called before the FFTW planner is created, then
#      FFTW's threads algorithms won't be registered or used at all.
#      (Previously, we called fftw_cleanup, but this invalidated existing
#       plans, causing Base Julia issue #19892.)
function __init__()
    check_deps()
    stat = ccall((:fftw_init_threads, libfftw3), Int32, ())
    statf = ccall((:fftwf_init_threads, libfftw3f), Int32, ())
    if stat == 0 || statf == 0
        error("could not initialize FFTW threads")
    end
end

include("fft.jl")
include("dct.jl")

end # module
