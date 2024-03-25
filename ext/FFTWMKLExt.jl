module FFTWMKLExt

using FFTW

# If we're using MKL, load it in and set library paths appropriately.
@static if FFTW.fftw_provider == "mkl"
    import MKL_jll
    FFTW.libfftw3[] = MKL_jll.libmkl_rt_path
    FFTW.libfftw3f[] = MKL_jll.libmkl_rt_path
    FFTW.version[] = VersionNumber(split(unsafe_string(cglobal(
                                    (:fftw_version,FFTW.libfftw3[]), UInt8)), ['-', ' '])[2])
end

function __init__()
    # Hook FFTW threads up to our partr runtime, and re-assign the
    # libfftw3{,f} refs at runtime, since we may have relocated and
    # changed the path to the library since the last time we precompiled.
    @static if FFTW.fftw_provider == "mkl"
        FFTW.libfftw3[] = MKL_jll.libmkl_rt_path
        FFTW.libfftw3f[] = MKL_jll.libmkl_rt_path
        FFTW.version[] = VersionNumber(split(unsafe_string(cglobal(
                                (:fftw_version,FFTW.libfftw3[]), UInt8)), ['-', ' '])[2])
    end
end

end
