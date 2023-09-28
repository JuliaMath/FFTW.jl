module FFTW

using LinearAlgebra, Reexport, Preferences
@reexport using AbstractFFTs
using Base.Threads

import AbstractFFTs: Plan, ScaledPlan,
                     fft, ifft, bfft, fft!, ifft!, bfft!,
                     plan_fft, plan_ifft, plan_bfft, plan_fft!, plan_ifft!, plan_bfft!,
                     rfft, irfft, brfft, plan_rfft, plan_irfft, plan_brfft,
                     fftshift, ifftshift,
                     rfft_output_size, brfft_output_size,
                     plan_inv, normalization

export dct, idct, dct!, idct!, plan_dct, plan_idct, plan_dct!, plan_idct!

include("providers.jl")

function __init__()
    # If someone is trying to set the provider via the old environment variable, warn them that they
    # should instead use `set_provider!()` instead.
    if haskey(ENV, "JULIA_FFTW_PROVIDER")
        Base.depwarn("JULIA_FFTW_PROVIDER is deprecated; use FFTW.set_provider!() instead", :JULIA_FFTW_PROVIDER)
    end

    # Hook FFTW threads up to our partr runtime, and re-assign the
    # libfftw3{,f} refs at runtime, since we may have relocated and
    # changed the path to the library since the last time we precompiled.
    @static if fftw_provider == "fftw"
        libfftw3[] = FFTW_jll.libfftw3_path
        libfftw3f[] = FFTW_jll.libfftw3f_path
        fftw_init_threads()
    end
    @static if fftw_provider == "mkl"
        libfftw3[] = MKL_jll.libmkl_rt_path
        libfftw3f[] = MKL_jll.libmkl_rt_path
    end
end

# most FFTW calls other than fftw_execute should be protected by a lock to be thread-safe
const fftwlock = ReentrantLock()

# protect thread-safety of expressions or whole functions by fftwlock:
macro exclusive(ex)
    if Meta.isexpr(ex, :function) || (Meta.isexpr(ex, :(=)) && Meta.isexpr(ex.args[1], :call))
        newbody = quote
            lock(fftwlock)
            try
                $(esc(ex.args[2]))
            finally
                unlock(fftwlock)
                destroy_deferred() # deallocate plans queued to be destroyed while we held the lock
            end
        end
        Expr(:function, esc(ex.args[1]), newbody)
    else
        return quote
            lock(fftwlock)
            try
                $(esc(ex))
            finally
                unlock(fftwlock)
                destroy_deferred() # deallocate plans queued to be destroyed while we held the lock
            end
        end
    end
end

include("fft.jl")
include("dct.jl")

include("precompile.jl")
_precompile_()

end # module
