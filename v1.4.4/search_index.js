var documenterSearchIndex = {"docs":
[{"location":"#FFTW.jl","page":"Home","title":"FFTW.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides bindings to the FFTW library for fast Fourier transforms.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is available for Julia versions 1.0 and later. To install it, run","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"FFTW\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"from the Julia REPL.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Users with a build of Julia based on Intel's Math Kernel Library (MKL) can use MKL for FFTs by setting a preference in their top-level project by either using the FFTW.set_provider!() method, or by directly setting the preference using Preferences.jl.  Note that this choice will be recorded for the current project, and other projects that wish to use MKL for FFTs should also set that same preference. Note further that MKL provides only a subset of the functionality provided by FFTW. See Intel's documentation for more information about potential differences or gaps in functionality.","category":"page"},{"location":"fft/#Fourier-Transforms","page":"API","title":"Fourier Transforms","text":"","category":"section"},{"location":"fft/","page":"API","title":"API","text":"This package extends the functionality provided by AbstractFFTs. To learn more about those functions, consult that package's documentation.","category":"page"},{"location":"fft/","page":"API","title":"API","text":"The following functions are unique to this package.","category":"page"},{"location":"fft/","page":"API","title":"API","text":"FFTW.r2r\nFFTW.r2r!\nFFTW.plan_r2r\nFFTW.plan_r2r!\nFFTW.dct\nFFTW.idct\nFFTW.dct!\nFFTW.idct!\nFFTW.plan_dct\nFFTW.plan_idct\nFFTW.plan_dct!\nFFTW.plan_idct!","category":"page"},{"location":"fft/#FFTW.r2r","page":"API","title":"FFTW.r2r","text":"r2r(A, kind [, dims])\n\nPerforms a multidimensional real-input/real-output (r2r) transform of type kind of the array A, as defined in the FFTW manual. kind specifies either a discrete cosine transform of various types (FFTW.REDFT00, FFTW.REDFT01, FFTW.REDFT10, or FFTW.REDFT11), a discrete sine transform of various types (FFTW.RODFT00, FFTW.RODFT01, FFTW.RODFT10, or FFTW.RODFT11), a real-input DFT with halfcomplex-format output (FFTW.R2HC and its inverse FFTW.HC2R), or a discrete Hartley transform (FFTW.DHT).  The kind argument may be an array or tuple in order to specify different transform types along the different dimensions of A; kind[end] is used for any unspecified dimensions.  See the FFTW manual for precise definitions of these transform types, at http://www.fftw.org/doc.\n\nThe optional dims argument specifies an iterable subset of dimensions (e.g. an integer, range, tuple, or array) to transform along. kind[i] is then the transform type for dims[i], with kind[end] being used for i > length(kind).\n\nSee also plan_r2r to pre-plan optimized r2r transforms.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.r2r!","page":"API","title":"FFTW.r2r!","text":"r2r!(A, kind [, dims])\n\nSame as r2r, but operates in-place on A, which must be an array of real or complex floating-point numbers.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.plan_r2r","page":"API","title":"FFTW.plan_r2r","text":"plan_r2r(A, kind [, dims [, flags [, timelimit]]])\n\nPre-plan an optimized r2r transform, similar to plan_fft except that the transforms (and the first three arguments) correspond to r2r and r2r!, respectively.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.plan_r2r!","page":"API","title":"FFTW.plan_r2r!","text":"plan_r2r!(A, kind [, dims [, flags [, timelimit]]])\n\nSimilar to plan_fft, but corresponds to r2r!.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.dct","page":"API","title":"FFTW.dct","text":"dct(A [, dims])\n\nPerforms a multidimensional type-II discrete cosine transform (DCT) of the array A, using the unitary normalization of the DCT. The optional dims argument specifies an iterable subset of dimensions (e.g. an integer, range, tuple, or array) to transform along.  Most efficient if the size of A along the transformed dimensions is a product of small primes; see nextprod. See also plan_dct for even greater efficiency.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.idct","page":"API","title":"FFTW.idct","text":"idct(A [, dims])\n\nComputes the multidimensional inverse discrete cosine transform (DCT) of the array A (technically, a type-III DCT with the unitary normalization). The optional dims argument specifies an iterable subset of dimensions (e.g. an integer, range, tuple, or array) to transform along.  Most efficient if the size of A along the transformed dimensions is a product of small primes; see nextprod.  See also plan_idct for even greater efficiency.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.dct!","page":"API","title":"FFTW.dct!","text":"dct!(A [, dims])\n\nSame as dct, except that it operates in-place on A, which must be an array of real or complex floating-point values.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.idct!","page":"API","title":"FFTW.idct!","text":"idct!(A [, dims])\n\nSame as idct, but operates in-place on A.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.plan_dct","page":"API","title":"FFTW.plan_dct","text":"plan_dct(A [, dims [, flags [, timelimit]]])\n\nPre-plan an optimized discrete cosine transform (DCT), similar to plan_fft except producing a function that computes dct. The first two arguments have the same meaning as for dct.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.plan_idct","page":"API","title":"FFTW.plan_idct","text":"plan_idct(A [, dims [, flags [, timelimit]]])\n\nPre-plan an optimized inverse discrete cosine transform (DCT), similar to plan_fft except producing a function that computes idct. The first two arguments have the same meaning as for idct.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.plan_dct!","page":"API","title":"FFTW.plan_dct!","text":"plan_dct!(A [, dims [, flags [, timelimit]]])\n\nSame as plan_dct, but operates in-place on A.\n\n\n\n\n\n","category":"function"},{"location":"fft/#FFTW.plan_idct!","page":"API","title":"FFTW.plan_idct!","text":"plan_idct!(A [, dims [, flags [, timelimit]]])\n\nSame as plan_idct, but operates in-place on A.\n\n\n\n\n\n","category":"function"}]
}
