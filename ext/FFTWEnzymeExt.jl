module FFTWEnzymeExt

using Enzyme, FFTW

Enzyme.EnzymeRules.inactive_noinl(typeof(assert_applicable), x...) = true

end # module

