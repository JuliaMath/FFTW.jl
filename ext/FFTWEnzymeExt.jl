module FFTWEnzymeExt

using Enzyme, FFTW

Enzyme.EnzymeRules.inactive_noinl(::typeof(FFTW.assert_applicable), x...) = true

Enzyme.EnzymeRules.inactive_noinl(::typeof(FFTW.unsafe_set_timelimit), x...) = true

end # module

