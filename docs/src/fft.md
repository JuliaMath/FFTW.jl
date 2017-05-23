# Fourier Transforms

```@docs
FFTW.fft
FFTW.fft!
FFTW.ifft
FFTW.ifft!
FFTW.bfft
FFTW.bfft!
FFTW.plan_fft
FFTW.plan_ifft
FFTW.plan_bfft
FFTW.plan_fft!
FFTW.plan_ifft!
FFTW.plan_bfft!
FFTW.rfft
FFTW.irfft
FFTW.brfft
FFTW.plan_rfft
FFTW.plan_brfft
FFTW.plan_irfft
FFTW.dct
FFTW.dct!
FFTW.idct
FFTW.idct!
FFTW.plan_dct
FFTW.plan_dct!
FFTW.plan_idct
FFTW.plan_idct!
FFTW.fftshift(::Any)
FFTW.fftshift(::Any, ::Any)
FFTW.ifftshift
```

The following functions are not exported from the package and thus must be qualified
with the `FFTW.` prefix on use.

```@docs
FFTW.r2r
FFTW.r2r!
FFTW.plan_r2r
FFTW.plan_r2r!
```
