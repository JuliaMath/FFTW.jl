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

## Note

These functions were formerly part of Base Julia.
Should any name conflicts occur on Julia versions which include these functions,
try adding

```julia
importall FFTW
```

to the top of your file.
If the problem persists, please open an issue on this package's repository.
