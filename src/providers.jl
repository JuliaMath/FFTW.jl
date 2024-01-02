const valid_fftw_providers = let
    # Hardcoded list of supported platforms
    # In principle, we could check FFTW_jll.is_available() and MKL_jll.is_available()
    # but then we would have to load MKL_jll which we want to avoid (lazy artifacts!)
    platforms_providers = Dict(
        Base.BinaryPlatforms.Platform("aarch64", "macos") => ("fftw",),
        Base.BinaryPlatforms.Platform("aarch64", "linux"; libc = "glibc") => ("fftw",),
        Base.BinaryPlatforms.Platform("aarch64", "linux"; libc = "musl") => ("fftw",),
        Base.BinaryPlatforms.Platform("armv6l", "linux"; libc = "glibc", call_abi = "eabihf") => ("fftw",),
        Base.BinaryPlatforms.Platform("armv6l", "linux"; libc = "musl", call_abi = "eabihf") => ("fftw",),
        Base.BinaryPlatforms.Platform("armv7l", "linux"; libc = "glibc", call_abi = "eabihf") => ("fftw",),
        Base.BinaryPlatforms.Platform("armv7l", "linux"; libc = "musl", call_abi = "eabihf") => ("fftw",),
        Base.BinaryPlatforms.Platform("i686", "linux"; libc = "glibc") => ("fftw", "mkl"),
        Base.BinaryPlatforms.Platform("i686", "linux"; libc = "musl") => ("fftw",),
        Base.BinaryPlatforms.Platform("i686", "windows") => ("fftw", "mkl"),
        Base.BinaryPlatforms.Platform("powerpc64le", "linux"; libc = "glibc") => ("fftw",),
        Base.BinaryPlatforms.Platform("x86_64", "macos") => ("fftw", "mkl"),
        Base.BinaryPlatforms.Platform("x86_64", "linux"; libc = "glibc") => ("fftw",),
        Base.BinaryPlatforms.Platform("x86_64", "linux"; libc = "musl") => ("fftw",),
        Base.BinaryPlatforms.Platform("x86_64", "freebsd") => ("fftw",),
        Base.BinaryPlatforms.Platform("x86_64", "windows") => ("fftw", "mkl"),
    )
    Base.BinaryPlatforms.select_platform(platforms_providers, Base.BinaryPlatforms.HostPlatform())
end
if valid_fftw_providers === nothing
    error("no valid FFTW library available")
end

function get_provider()
    # Note: we CANNOT do something like have the `default` value be `get(ENV, "JULIA_FFTW_PROVIDER", "fftw")` here.
    # This is because the only way the Julia knows that a default has changed is if the values on-disk change; so
    # if your "default" value can be changed from the outside, you quickly run into cache invalidation issues.
    # So the default here _must_ be a constant.
    default_provider = first(valid_fftw_providers)

    # Load the preference
    provider = @load_preference("provider", default_provider)

    # Ensure the provider matches one of the ones we support
    if provider ∉ valid_fftw_providers
        @error("Invalid provider setting \"$(provider)\"; valid settings include [$(join(map(x -> '"' * x * '"', valid_fftw_providers), ", "))]")
        provider = default_provider
    end
    return provider
end

# Read in preferences, see if any users have requested a particular backend
const fftw_provider = get_provider()

# We'll initialize `libfftw3` here (in the conditionals below), and
# it will get overwritten again in `__init__()`.  This allows us to
# `ccall` at build time, and also be relocatable for PackageCompiler.
const libfftw3 = Ref{String}()
const libfftw3f = Ref{String}()

"""
    set_provider!(provider; export_prefs::Bool = false)

Convenience wrapper for setting the FFT provider.  Valid values include `"fftw"`, `"mkl"`.
Also supports `Preferences` sentinel values `nothing` and `missing`; see the docstring for
`Preferences.set_preferences!()` for more information on what these values mean.
"""
function set_provider!(provider; export_prefs::Bool = false)
    if provider !== nothing && provider !== missing && provider ∉ valid_fftw_providers
        throw(ArgumentError("Invalid provider value \"$(provider)\"; valid settings include [$(join(map(x -> '"' * x * '"', valid_fftw_providers), ", "))]"))
    end
    set_preferences!(@__MODULE__, "provider" => provider; export_prefs, force = true)
    if provider != fftw_provider
        # Re-fetch to get default values in the event that `nothing` or `missing` was passed in.
        provider = get_provider()
        @info("FFTW provider changed; restart Julia for this change to take effect", provider)
    end
end

# If we're using fftw_jll, load it in
@static if fftw_provider == "fftw"
    import FFTW_jll
    if !FFTW_jll.is_available()
        # more descriptive error message if FFTW is not available
        # (should not be possible to reach this)
        @error("FFTW library cannot be loaded: Run `FFTW.set_provider!(\"mkl\")` to switch to MKL")
    end
    libfftw3[] = FFTW_jll.libfftw3_path
    libfftw3f[] = FFTW_jll.libfftw3f_path

    # callback function that FFTW uses to launch `num` parallel
    # tasks (FFTW/fftw3#175):
    function spawnloop(f::Ptr{Cvoid}, fdata::Ptr{Cvoid}, elsize::Csize_t, num::Cint, callback_data::Ptr{Cvoid})
        @sync for i = 0:num-1
            Threads.@spawn ccall(f, Ptr{Cvoid}, (Ptr{Cvoid},), fdata + elsize*i)
        end
    end

    # If FFTW was built with threads, then they must be initialized before any FFTW planning routine.
    #   -- This initializes FFTW's threads support (defaulting to 1 thread).
    #      If this isn't called before the FFTW planner is created, then
    #      FFTW's threads algorithms won't be registered or used at all.
    #      (Previously, we called fftw_cleanup, but this invalidated existing
    #       plans, causing Base Julia issue #19892.)
    function fftw_init_threads()
        stat =  ccall((:fftw_init_threads,   libfftw3[]), Int32, ())
        statf = ccall((:fftwf_init_threads, libfftw3f[]), Int32, ())
        if stat == 0 || statf == 0
            error("could not initialize FFTW threads")
        end

        if nthreads() > 1
            cspawnloop = @cfunction(spawnloop, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t, Cint, Ptr{Cvoid}))
            ccall((:fftw_threads_set_callback,  libfftw3[]),  Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), cspawnloop, C_NULL)
            ccall((:fftwf_threads_set_callback, libfftw3f[]), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), cspawnloop, C_NULL)
        end
    end
end

# If we're using MKL, load it in and set library paths appropriately.
@static if fftw_provider == "mkl"
    import MKL_jll
    if !MKL_jll.is_available()
        # more descriptive error message if MKL is not available
        # (should not be possible to reach this)
        @error("MKL cannot be loaded: Run `FFTW.set_provider!(\"fftw\")` to switch to the FFTW library")
    end
    libfftw3[] = MKL_jll.libmkl_rt_path
    libfftw3f[] = MKL_jll.libmkl_rt_path
    const _last_num_threads = Ref(Cint(1))
end
