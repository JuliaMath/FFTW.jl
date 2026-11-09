# Fourier Transforms

This package extends the functionality provided by
[AbstractFFTs](https://github.com/JuliaMath/AbstractFFTs.jl).
To learn more about those functions, consult that package's
[documentation](https://juliamath.github.io/AbstractFFTs.jl/stable/api/).

## Planning with a thread-count override

FFTW's implementations of `plan_fft`, `plan_fft!`, `plan_bfft`, `plan_bfft!`,
`plan_rfft` and `plan_brfft` accept the keyword `num_threads`. This is an
FFTW-specific extension to the generic AbstractFFTs planning interface.

For example, create and apply an in-place complex FFT plan using a thread count
of two:

```jldoctest
julia> using FFTW

julia> previous_threads = FFTW.get_num_threads();

julia> x = ComplexF64[1, 2, 3, 4];

julia> p = plan_fft!(x; flags=FFTW.ESTIMATE, num_threads=2);

julia> FFTW.get_num_threads() == previous_threads
true

julia> p * x ≈ ComplexF64[10, -2+2im, -2, -2-2im]
true
```

`num_threads` overrides the planner thread setting while this plan is created;
the previous setting is restored afterwards. Omitting it (or passing `nothing`)
uses the current setting from `FFTW.get_num_threads()`, which can be changed
with `FFTW.set_num_threads(n)` for subsequent plans. The keyword does not change
existing plans or Julia's `Threads.nthreads()`. It is a planning choice, not a
guarantee of a speedup for a particular transform. `plan_fft!` creates an
in-place plan, so applying `p` overwrites `x` in the example.

## FFTW-specific functions

The following functions are unique to this package.

```@docs
FFTW.r2r
FFTW.r2r!
FFTW.plan_r2r
FFTW.plan_r2r!
FFTW.dct
FFTW.idct
FFTW.dct!
FFTW.idct!
FFTW.plan_dct
FFTW.plan_idct
FFTW.plan_dct!
FFTW.plan_idct!
FFTW.set_provider!
```
