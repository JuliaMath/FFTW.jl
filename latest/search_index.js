var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#FFTW.jl-1",
    "page": "Home",
    "title": "FFTW.jl",
    "category": "section",
    "text": "This package provides bindings to the FFTW library for fast Fourier transforms."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is available for Julia versions 0.6 and up. To install it, runPkg.add(\"FFTW\")from the Julia REPL."
},

{
    "location": "index.html#Note-1",
    "page": "Home",
    "title": "Note",
    "category": "section",
    "text": "These functions were formerly part of Base Julia. Should any name conflicts occur on Julia versions which include these functions, try addingimportall FFTWto the top of your file. If the problem persists, please open an issue on this package's repository."
},

{
    "location": "fft.html#",
    "page": "Fourier Transforms",
    "title": "Fourier Transforms",
    "category": "page",
    "text": ""
},

{
    "location": "fft.html#FFTW.fft",
    "page": "Fourier Transforms",
    "title": "FFTW.fft",
    "category": "Function",
    "text": "fft(A [, dims])\n\nPerforms a multidimensional FFT of the array A. The optional dims argument specifies an iterable subset of dimensions (e.g. an integer, range, tuple, or array) to transform along. Most efficient if the size of A along the transformed dimensions is a product of small primes; see nextprod(). See also plan_fft() for even greater efficiency.\n\nA one-dimensional FFT computes the one-dimensional discrete Fourier transform (DFT) as defined by\n\noperatornameDFT(A)k =\n  sum_n=1^operatornamelength(A)\n  expleft(-ifrac2pi\n  (n-1)(k-1)operatornamelength(A) right) An\n\nA multidimensional FFT simply performs this operation along each transformed dimension of A.\n\nnote: Note\nJulia starts FFTW up with 1 thread by default. Higher performance is usually possible by increasing number of threads. Use FFTW.set_num_threads(Sys.CPU_CORES) to use as many threads as cores on your system.\nThis performs a multidimensional FFT by default. FFT libraries in other languages such as Python and Octave perform a one-dimensional FFT along the first non-singleton dimension of the array. This is worth noting while performing comparisons. For more details, refer to the Noteworthy Differences from other Languages section of the manual.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.fft!",
    "page": "Fourier Transforms",
    "title": "FFTW.fft!",
    "category": "Function",
    "text": "fft!(A [, dims])\n\nSame as fft, but operates in-place on A, which must be an array of complex floating-point numbers.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.ifft",
    "page": "Fourier Transforms",
    "title": "FFTW.ifft",
    "category": "Function",
    "text": "ifft(A [, dims])\n\nMultidimensional inverse FFT.\n\nA one-dimensional inverse FFT computes\n\noperatornameIDFT(A)k = frac1operatornamelength(A)\nsum_n=1^operatornamelength(A) expleft(+ifrac2pi (n-1)(k-1)\noperatornamelength(A) right) An\n\nA multidimensional inverse FFT simply performs this operation along each transformed dimension of A.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.ifft!",
    "page": "Fourier Transforms",
    "title": "FFTW.ifft!",
    "category": "Function",
    "text": "ifft!(A [, dims])\n\nSame as ifft, but operates in-place on A.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.bfft",
    "page": "Fourier Transforms",
    "title": "FFTW.bfft",
    "category": "Function",
    "text": "bfft(A [, dims])\n\nSimilar to ifft, but computes an unnormalized inverse (backward) transform, which must be divided by the product of the sizes of the transformed dimensions in order to obtain the inverse. (This is slightly more efficient than ifft because it omits a scaling step, which in some applications can be combined with other computational steps elsewhere.)\n\noperatornameBDFT(A)k = operatornamelength(A) operatornameIDFT(A)k\n\n\n\n"
},

{
    "location": "fft.html#FFTW.bfft!",
    "page": "Fourier Transforms",
    "title": "FFTW.bfft!",
    "category": "Function",
    "text": "bfft!(A [, dims])\n\nSame as bfft, but operates in-place on A.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_fft",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_fft",
    "category": "Function",
    "text": "plan_fft(A [, dims]; flags=FFTW.ESTIMATE;  timelimit=Inf)\n\nPre-plan an optimized FFT along given dimensions (dims) of arrays matching the shape and type of A.  (The first two arguments have the same meaning as for fft.) Returns an object P which represents the linear operator computed by the FFT, and which contains all of the information needed to compute fft(A, dims) quickly.\n\nTo apply P to an array A, use P * A; in general, the syntax for applying plans is much like that of matrices.  (A plan can only be applied to arrays of the same size as the A for which the plan was created.)  You can also apply a plan with a preallocated output array Â by calling A_mul_B!(Â, plan, A).  (For A_mul_B!, however, the input array A must be a complex floating-point array like the output Â.) You can compute the inverse-transform plan by inv(P) and apply the inverse plan with P \\ Â (the inverse plan is cached and reused for subsequent calls to inv or \\), and apply the inverse plan to a pre-allocated output array A with A_ldiv_B!(A, P, Â).\n\nThe flags argument is a bitwise-or of FFTW planner flags, defaulting to FFTW.ESTIMATE. e.g. passing FFTW.MEASURE or FFTW.PATIENT will instead spend several seconds (or more) benchmarking different possible FFT algorithms and picking the fastest one; see the FFTW manual for more information on planner flags.  The optional timelimit argument specifies a rough upper bound on the allowed planning time, in seconds. Passing FFTW.MEASURE or FFTW.PATIENT may cause the input array A to be overwritten with zeros during plan creation.\n\nplan_fft! is the same as plan_fft but creates a plan that operates in-place on its argument (which must be an array of complex floating-point numbers). plan_ifft and so on are similar but produce plans that perform the equivalent of the inverse transforms ifft and so on.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_ifft",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_ifft",
    "category": "Function",
    "text": "plan_ifft(A [, dims]; flags=FFTW.ESTIMATE;  timelimit=Inf)\n\nSame as plan_fft, but produces a plan that performs inverse transforms ifft.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_bfft",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_bfft",
    "category": "Function",
    "text": "plan_bfft(A [, dims]; flags=FFTW.ESTIMATE;  timelimit=Inf)\n\nSame as plan_fft, but produces a plan that performs an unnormalized backwards transform bfft.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_fft!",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_fft!",
    "category": "Function",
    "text": "plan_fft!(A [, dims]; flags=FFTW.ESTIMATE;  timelimit=Inf)\n\nSame as plan_fft, but operates in-place on A.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_ifft!",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_ifft!",
    "category": "Function",
    "text": "plan_ifft!(A [, dims]; flags=FFTW.ESTIMATE;  timelimit=Inf)\n\nSame as plan_ifft, but operates in-place on A.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_bfft!",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_bfft!",
    "category": "Function",
    "text": "plan_bfft!(A [, dims]; flags=FFTW.ESTIMATE;  timelimit=Inf)\n\nSame as plan_bfft, but operates in-place on A.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.rfft",
    "page": "Fourier Transforms",
    "title": "FFTW.rfft",
    "category": "Function",
    "text": "rfft(A [, dims])\n\nMultidimensional FFT of a real array A, exploiting the fact that the transform has conjugate symmetry in order to save roughly half the computational time and storage costs compared with fft. If A has size (n_1, ..., n_d), the result has size (div(n_1,2)+1, ..., n_d).\n\nThe optional dims argument specifies an iterable subset of one or more dimensions of A to transform, similar to fft. Instead of (roughly) halving the first dimension of A in the result, the dims[1] dimension is (roughly) halved in the same way.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.irfft",
    "page": "Fourier Transforms",
    "title": "FFTW.irfft",
    "category": "Function",
    "text": "irfft(A, d [, dims])\n\nInverse of rfft: for a complex array A, gives the corresponding real array whose FFT yields A in the first half. As for rfft, dims is an optional subset of dimensions to transform, defaulting to 1:ndims(A).\n\nd is the length of the transformed real array along the dims[1] dimension, which must satisfy div(d,2)+1 == size(A,dims[1]). (This parameter cannot be inferred from size(A) since both 2*size(A,dims[1])-2 as well as 2*size(A,dims[1])-1 are valid sizes for the transformed real array.)\n\n\n\n"
},

{
    "location": "fft.html#FFTW.brfft",
    "page": "Fourier Transforms",
    "title": "FFTW.brfft",
    "category": "Function",
    "text": "brfft(A, d [, dims])\n\nSimilar to irfft but computes an unnormalized inverse transform (similar to bfft), which must be divided by the product of the sizes of the transformed dimensions (of the real output array) in order to obtain the inverse transform.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_rfft",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_rfft",
    "category": "Function",
    "text": "plan_rfft(A [, dims]; flags=FFTW.ESTIMATE;  timelimit=Inf)\n\nPre-plan an optimized real-input FFT, similar to plan_fft except for rfft instead of fft. The first two arguments, and the size of the transformed result, are the same as for rfft.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_brfft",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_brfft",
    "category": "Function",
    "text": "plan_brfft(A, d [, dims]; flags=FFTW.ESTIMATE;  timelimit=Inf)\n\nPre-plan an optimized real-input unnormalized transform, similar to plan_rfft except for brfft instead of rfft. The first two arguments and the size of the transformed result, are the same as for brfft.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_irfft",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_irfft",
    "category": "Function",
    "text": "plan_irfft(A, d [, dims]; flags=FFTW.ESTIMATE;  timelimit=Inf)\n\nPre-plan an optimized inverse real-input FFT, similar to plan_rfft except for irfft and brfft, respectively. The first three arguments have the same meaning as for irfft.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.dct",
    "page": "Fourier Transforms",
    "title": "FFTW.dct",
    "category": "Function",
    "text": "dct(A [, dims])\n\nPerforms a multidimensional type-II discrete cosine transform (DCT) of the array A, using the unitary normalization of the DCT. The optional dims argument specifies an iterable subset of dimensions (e.g. an integer, range, tuple, or array) to transform along.  Most efficient if the size of A along the transformed dimensions is a product of small primes; see nextprod. See also plan_dct for even greater efficiency.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.dct!",
    "page": "Fourier Transforms",
    "title": "FFTW.dct!",
    "category": "Function",
    "text": "dct!(A [, dims])\n\nSame as dct!, except that it operates in-place on A, which must be an array of real or complex floating-point values.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.idct",
    "page": "Fourier Transforms",
    "title": "FFTW.idct",
    "category": "Function",
    "text": "idct(A [, dims])\n\nComputes the multidimensional inverse discrete cosine transform (DCT) of the array A (technically, a type-III DCT with the unitary normalization). The optional dims argument specifies an iterable subset of dimensions (e.g. an integer, range, tuple, or array) to transform along.  Most efficient if the size of A along the transformed dimensions is a product of small primes; see nextprod.  See also plan_idct for even greater efficiency.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.idct!",
    "page": "Fourier Transforms",
    "title": "FFTW.idct!",
    "category": "Function",
    "text": "idct!(A [, dims])\n\nSame as idct!, but operates in-place on A.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_dct",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_dct",
    "category": "Function",
    "text": "plan_dct(A [, dims [, flags [, timelimit]]])\n\nPre-plan an optimized discrete cosine transform (DCT), similar to plan_fft except producing a function that computes dct. The first two arguments have the same meaning as for dct.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_dct!",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_dct!",
    "category": "Function",
    "text": "plan_dct!(A [, dims [, flags [, timelimit]]])\n\nSame as plan_dct, but operates in-place on A.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_idct",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_idct",
    "category": "Function",
    "text": "plan_idct(A [, dims [, flags [, timelimit]]])\n\nPre-plan an optimized inverse discrete cosine transform (DCT), similar to plan_fft except producing a function that computes idct. The first two arguments have the same meaning as for idct.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_idct!",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_idct!",
    "category": "Function",
    "text": "plan_idct!(A [, dims [, flags [, timelimit]]])\n\nSame as plan_idct, but operates in-place on A.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.fftshift-Tuple{Any}",
    "page": "Fourier Transforms",
    "title": "FFTW.fftshift",
    "category": "Method",
    "text": "fftshift(x)\n\nSwap the first and second halves of each dimension of x.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.fftshift-Tuple{Any,Any}",
    "page": "Fourier Transforms",
    "title": "FFTW.fftshift",
    "category": "Method",
    "text": "fftshift(x,dim)\n\nSwap the first and second halves of the given dimension or iterable of dimensions of array x.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.ifftshift",
    "page": "Fourier Transforms",
    "title": "FFTW.ifftshift",
    "category": "Function",
    "text": "ifftshift(x, [dim])\n\nUndoes the effect of fftshift.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.r2r",
    "page": "Fourier Transforms",
    "title": "FFTW.r2r",
    "category": "Function",
    "text": "r2r(A, kind [, dims])\n\nPerforms a multidimensional real-input/real-output (r2r) transform of type kind of the array A, as defined in the FFTW manual. kind specifies either a discrete cosine transform of various types (FFTW.REDFT00, FFTW.REDFT01, FFTW.REDFT10, or FFTW.REDFT11), a discrete sine transform of various types (FFTW.RODFT00, FFTW.RODFT01, FFTW.RODFT10, or FFTW.RODFT11), a real-input DFT with halfcomplex-format output (FFTW.R2HC and its inverse FFTW.HC2R), or a discrete Hartley transform (FFTW.DHT).  The kind argument may be an array or tuple in order to specify different transform types along the different dimensions of A; kind[end] is used for any unspecified dimensions.  See the FFTW manual for precise definitions of these transform types, at http://www.fftw.org/doc.\n\nThe optional dims argument specifies an iterable subset of dimensions (e.g. an integer, range, tuple, or array) to transform along. kind[i] is then the transform type for dims[i], with kind[end] being used for i > length(kind).\n\nSee also plan_r2r to pre-plan optimized r2r transforms.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.r2r!",
    "page": "Fourier Transforms",
    "title": "FFTW.r2r!",
    "category": "Function",
    "text": "r2r!(A, kind [, dims])\n\nSame as r2r, but operates in-place on A, which must be an array of real or complex floating-point numbers.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_r2r",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_r2r",
    "category": "Function",
    "text": "plan_r2r(A, kind [, dims [, flags [, timelimit]]])\n\nPre-plan an optimized r2r transform, similar to plan_fft except that the transforms (and the first three arguments) correspond to r2r and r2r!, respectively.\n\n\n\n"
},

{
    "location": "fft.html#FFTW.plan_r2r!",
    "page": "Fourier Transforms",
    "title": "FFTW.plan_r2r!",
    "category": "Function",
    "text": "plan_r2r!(A, kind [, dims [, flags [, timelimit]]])\n\nSimilar to plan_fft, but corresponds to r2r!.\n\n\n\n"
},

{
    "location": "fft.html#Fourier-Transforms-1",
    "page": "Fourier Transforms",
    "title": "Fourier Transforms",
    "category": "section",
    "text": "FFTW.fft\nFFTW.fft!\nFFTW.ifft\nFFTW.ifft!\nFFTW.bfft\nFFTW.bfft!\nFFTW.plan_fft\nFFTW.plan_ifft\nFFTW.plan_bfft\nFFTW.plan_fft!\nFFTW.plan_ifft!\nFFTW.plan_bfft!\nFFTW.rfft\nFFTW.irfft\nFFTW.brfft\nFFTW.plan_rfft\nFFTW.plan_brfft\nFFTW.plan_irfft\nFFTW.dct\nFFTW.dct!\nFFTW.idct\nFFTW.idct!\nFFTW.plan_dct\nFFTW.plan_dct!\nFFTW.plan_idct\nFFTW.plan_idct!\nFFTW.fftshift(::Any)\nFFTW.fftshift(::Any, ::Any)\nFFTW.ifftshiftThe following functions are not exported from the package and thus must be qualified with the FFTW. prefix on use.FFTW.r2r\nFFTW.r2r!\nFFTW.plan_r2r\nFFTW.plan_r2r!"
},

{
    "location": "dsp.html#",
    "page": "Signal Processing",
    "title": "Signal Processing",
    "category": "page",
    "text": ""
},

{
    "location": "dsp.html#FFTW.filt",
    "page": "Signal Processing",
    "title": "FFTW.filt",
    "category": "Function",
    "text": "filt(b, a, x, [si])\n\nApply filter described by vectors a and b to vector x, with an optional initial filter state vector si (defaults to zeros).\n\n\n\n"
},

{
    "location": "dsp.html#FFTW.filt!",
    "page": "Signal Processing",
    "title": "FFTW.filt!",
    "category": "Function",
    "text": "filt!(out, b, a, x, [si])\n\nSame as filt but writes the result into the out argument, which may alias the input x to modify it in-place.\n\n\n\n"
},

{
    "location": "dsp.html#FFTW.deconv",
    "page": "Signal Processing",
    "title": "FFTW.deconv",
    "category": "Function",
    "text": "deconv(b,a) -> c\n\nConstruct vector c such that b = conv(a,c) + r. Equivalent to polynomial division.\n\n\n\n"
},

{
    "location": "dsp.html#FFTW.conv",
    "page": "Signal Processing",
    "title": "FFTW.conv",
    "category": "Function",
    "text": "conv(u,v)\n\nConvolution of two vectors. Uses FFT algorithm.\n\n\n\n"
},

{
    "location": "dsp.html#FFTW.conv2",
    "page": "Signal Processing",
    "title": "FFTW.conv2",
    "category": "Function",
    "text": "conv2(u,v,A)\n\n2-D convolution of the matrix A with the 2-D separable kernel generated by the vectors u and v. Uses 2-D FFT algorithm.\n\n\n\nconv2(B,A)\n\n2-D convolution of the matrix B with the matrix A. Uses 2-D FFT algorithm.\n\n\n\n"
},

{
    "location": "dsp.html#FFTW.xcorr",
    "page": "Signal Processing",
    "title": "FFTW.xcorr",
    "category": "Function",
    "text": "xcorr(u,v)\n\nCompute the cross-correlation of two vectors.\n\n\n\n"
},

{
    "location": "dsp.html#Signal-Processing-1",
    "page": "Signal Processing",
    "title": "Signal Processing",
    "category": "section",
    "text": "It is expected that these functions will be moved to a different package.FFTW.filt\nFFTW.filt!\nFFTW.deconv\nFFTW.conv\nFFTW.conv2\nFFTW.xcorr"
},

]}
