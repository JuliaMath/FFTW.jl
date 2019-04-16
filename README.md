# FFTW.jl

[![Travis](https://travis-ci.org/JuliaMath/FFTW.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/FFTW.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/hofbdbyt287qn49s/branch/master?svg=true)](https://ci.appveyor.com/project/ararslan/fftw-jl/branch/master)
[![Coveralls](https://coveralls.io/repos/github/JuliaMath/FFTW.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/FFTW.jl?branch=master)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/FFTW.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaMath.github.io/FFTW.jl/latest)

This package provides Julia bindings to the [FFTW](http://www.fftw.org/) library for
fast Fourier transforms (FFTs), as well as functionality useful for signal processing.
These functions were formerly a part of Base Julia.

## Usage and documentation

```julia
]add FFTW
using FFTW
fft([0; 1; 2; 1])
```
returns
```
4-element Array{Complex{Float64},1}:
  4.0 + 0.0im
 -2.0 + 0.0im
  0.0 + 0.0im
 -2.0 + 0.0im
```

The documentation of generic FFT functionality can be found in the [AbstractFFTs.jl package](https://juliamath.github.io/AbstractFFTs.jl/stable/api/#Public-Interface-1). Additional functionalities supported by the FFTW library are documented in the present package.

## MKL

Alternatively, the FFTs in Intel's Math Kernel Library (MKL) can be used
by setting an environment variable `JULIA_FFTW_PROVIDER` to `MKL` and running
`Pkg.build("FFTW")`.   If Julia was built with MKL, then Julia's MKL will be
used for FFTs; otherwise the [Conda.jl package](https://github.com/JuliaPy/Conda.jl)
will be used to download MKL via the [mkl_fft Anaconda package](https://github.com/IntelPython/mkl_fft).
Setting this environment variable only needs to be done for the first build of the package;
after that, the package will remember to use MKL when building and updating.
Note however that MKL provides only a subset of the functionality provided by FFTW. See
Intel's [documentation](https://software.intel.com/en-us/mkl-developer-reference-c-using-fftw3-wrappers)
for more information about potential differences or gaps in functionality.


## License

The FFTW library will be downloaded on versions of Julia where it is no longer distributed
as part of Julia.
Note that FFTW is licensed under GPLv2 or higher (see
[its license file](http://www.fftw.org/doc/License-and-Copyright.html)), but the bindings
to the library in this package, FFTW.jl, are licensed under MIT.
This means that code using the FFTW library via the FFTW.jl bindings is subject to FFTW's
licensing terms.
Code using alternative implementations of the FFTW API, such as
[MKL's FFTW3 interface](https://software.intel.com/en-us/mkl-developer-reference-c-fftw3-interface-to-intel-math-kernel-library)
are instead subject to the alternative's license.
If you distribute a derived or combined work, i.e. a program that links to and is distributed
with the FFTW library, then that distribution falls under the terms of the GPL.
If you just distribute source code that links to FFTW.jl, and users have to download FFTW
or MKL to provide the backend, then the GPL probably doesn't have much effect on you.
