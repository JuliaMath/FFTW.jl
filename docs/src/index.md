# FFTW.jl

This package provides bindings to the [FFTW](http://www.fftw.org/) library for
fast Fourier transforms.

## Installation

The package is available for Julia versions 1.0 and later.
To install it, run

```julia
using Pkg
Pkg.add("FFTW")
```

from the Julia REPL.

Users with a build of Julia based on Intel's Math Kernel Library (MKL) can use MKL for FFTs by setting a preference in
their top-level project by either using the [`FFTW.set_provider!()`](@ref) method, or by directly setting the preference
using [`Preferences.jl`](https://github.com/JuliaPackaging/Preferences.jl).  Note that this choice will be recorded for
the current project, and other projects that wish to use MKL for FFTs should also set that same preference. Note further
that MKL provides only a subset of the functionality provided by FFTW. See Intel's
[documentation](https://software.intel.com/en-us/mkl-developer-reference-c-using-fftw3-wrappers) for more information
about potential differences or gaps in functionality.
