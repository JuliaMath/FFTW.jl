# FFTW.jl

This package provides Julia bindings to the [FFTW](http://www.fftw.org/) library for
fast Fourier transforms, as well as functionality useful for signal processing.
These functions were formerly a part of Base Julia.

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

[![Travis](https://travis-ci.org/JuliaMath/FFTW.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/FFTW.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/hofbdbyt287qn49s/branch/master?svg=true)](https://ci.appveyor.com/project/ararslan/fftw-jl/branch/master)
[![Coveralls](https://coveralls.io/repos/github/JuliaMath/FFTW.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/FFTW.jl?branch=master)

Documentation:
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/FFTW.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaMath.github.io/FFTW.jl/latest)
