# FFTW.jl

This package provides bindings to the [FFTW](http://www.fftw.org/) library for
fast Fourier transforms.

## Installation

The package is available for Julia versions 0.6 and up.
To install it, run

```julia
Pkg.add("FFTW")
```

from the Julia REPL.

Users with a build of Julia based on Intel's Math Kernel Library (MKL) can take use MKL
for FFTs by setting an environment variable `JULIA_FFTW_PROVIDER` to `MKL` and running
`Pkg.build("FFTW")`.
Setting this environment variable only needs to be done for the first build of the package;
after that, the package will remember to use MKL when building and updating.
Note however that MKL provides only a subset of the functionality provided by FFTW. See
Intel's [documentation](https://software.intel.com/en-us/mkl-developer-reference-c-using-fftw3-wrappers)
for more information about potential differences or gaps in functionality.

## Note

These functions were formerly part of Base Julia.
Should any name conflicts occur on Julia versions which include these functions,
try adding

```julia
importall FFTW
```

to the top of your file.
If the problem persists, please open an issue on this package's repository.
