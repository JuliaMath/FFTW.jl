module FFTW

using LinearAlgebra, Reexport
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

# Use Julia partr threading backend if present
@static if fftw_vendor == :fftw && isdefined(Threads, Symbol("@spawn"))
    # callback function that FFTW uses to launch `num` parallel
    # tasks (FFTW/fftw3#175):
    function spawnloop(f::Ptr{Cvoid}, fdata::Ptr{Cvoid}, elsize::Csize_t, num::Cint, callback_data::Ptr{Cvoid})
        @sync for i = 0:num-1
            Threads.@spawn ccall(f, Ptr{Cvoid}, (Ptr{Cvoid},), fdata + elsize*i)
        end
    end
end

# If FFTW was built with threads, then they must be initialized before any FFTW planning routine.
#   -- This initializes FFTW's threads support (defaulting to 1 thread).
#      If this isn't called before the FFTW planner is created, then
#      FFTW's threads algorithms won't be registered or used at all.
#      (Previously, we called fftw_cleanup, but this invalidated existing
#       plans, causing Base Julia issue #19892.)
function __init__()
    check_deps()
    stat =  ccall((:fftw_init_threads,   libfftw3), Int32, ())
    statf = ccall((:fftwf_init_threads, libfftw3f), Int32, ())
    if stat == 0 || statf == 0
        error("could not initialize FFTW threads")
    end
    @static if fftw_vendor == :fftw
        if Threads.nthreads() > 1 # number of Julia threads is set when Julia is launched
            ccall((:fftw_make_planner_thread_safe,  libfftw3),  Cvoid, ())
            ccall((:fftwf_make_planner_thread_safe, libfftw3f), Cvoid, ())
        end
        @static if isdefined(Threads, Symbol("@spawn"))
            if Threads.nthreads() > 1 # partr will give us our threads
                cspawnloop = @cfunction(spawnloop, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t, Cint, Ptr{Cvoid}))
                ccall((:fftw_threads_set_callback,  libfftw3),  Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), cspawnloop, C_NULL)
                ccall((:fftwf_threads_set_callback, libfftw3f), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), cspawnloop, C_NULL)
                set_num_threads(Threads.nthreads() * 4) # spawn more tasks than threads to help load-balancing
            end
        end
    end
end

include("fft.jl")
include("dct.jl")

end # module
