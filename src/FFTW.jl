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

function fftw_init_check()
    # If someone is trying to set the provider via the old environment variable, warn them that they
    # should instead use `set_provider!()` instead.
    if haskey(ENV, "JULIA_FFTW_PROVIDER")
        Base.depwarn("JULIA_FFTW_PROVIDER is deprecated; use FFTW.set_provider!() instead", :JULIA_FFTW_PROVIDER)
    end

    # Hook FFTW threads up to our partr runtime
    @static if fftw_provider == "fftw"
        fftw_init_threads()
    end
end

if VERSION >= v"1.11.0"
# This can be be deleted once FFTW_jll is upgraded to the real lazy jll code, to get the real benefits of this mess
mutable struct FakeLazyLibrary
    reallibrary::Symbol
    on_load_callback
    @atomic h::Ptr{Cvoid}
end
import Libdl: LazyLibrary, dlopen
function dlopen(lib::FakeLazyLibrary)
    h = @atomic :monotonic lib.h
    h != C_NULL && return h
    @lock fftwlock begin
        h = @atomic :monotonic lib.h
        h != C_NULL && return h
        h = dlopen(getglobal(FFTW, lib.reallibrary))
        lib.on_load_callback()
        @atomic :release lib.h = h
    end
    return h
end

@static if fftw_provider == "fftw"
    import FFTW_jll: libfftw3 as libfftw3_no_init,
                     libfftw3f as libfftw3f_no_init
elseif fftw_provider == "mkl"
    import MKL_jll: libmkl_rt as libfftw3_no_init,
                    libmkl_rt as libfftw3f_no_init
end
const libfftw3 = FakeLazyLibrary(:libfftw3_no_init, fftw_init_check, C_NULL)
const libfftw3f = FakeLazyLibrary(:libfftw3f_no_init, fftw_init_check, C_NULL)

else
@static if fftw_provider == "fftw"
    import FFTW_jll: libfftw3_path as libfftw3_no_init,
                     libfftw3f_path as libfftw3f_no_init,
                     libfftw3_path as libfftw3,
                     libfftw3f_path as libfftw3f
elseif fftw_provider == "mkl"
    import MKL_jll: libmkl_rt_path as libfftw3_no_init,
                    libmkl_rt_path as libfftw3f_no_init,
                    libmkl_rt_path as libfftw3,
                    libmkl_rt_path as libfftw3f
end
function __init__()
    fftw_init_check()
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
