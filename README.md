# FFTW.jl

This package provides Julia bindings to the [FFTW](http://www.fftw.org/) library for
fast Fourier transforms, as well as functionality useful for signal processing.
These functions were formerly a part of Base Julia.

Users with a build of Julia based on Intel's Math Kernel Library (MKL) can take use MKL
for FFTs by setting an environment variable `JULIA_FFTW_PROVIDER` to `MKL` and running
`Pkg.build("FFTW")`.
Setting this environment variable only needs to be done for the first build of the package;
after that, the package will remember to use MKL when building and updating.
Note however that MKL provides only a subset of the functionality provided by FFTW. See
Intel's [documentation](https://software.intel.com/en-us/mkl-developer-reference-c-using-fftw3-wrappers)
for more information about potential differences or gaps in functionality.

The FFTW library will be downloaded on versions of Julia where it is no longer distributed
as part of Julia.
Note that FFTW is licensed under GPLv2 or higher (see
[its license file](http://www.fftw.org/doc/License-and-Copyright.html)), but the bindings
here are licensed under MIT.

[![Travis](https://travis-ci.org/JuliaMath/FFTW.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/FFTW.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/hofbdbyt287qn49s/branch/master?svg=true)](https://ci.appveyor.com/project/ararslan/fftw-jl/branch/master)
[![Coveralls](https://coveralls.io/repos/github/JuliaMath/FFTW.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/FFTW.jl?branch=master)

Documentation:
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/FFTW.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaMath.github.io/FFTW.jl/latest)
