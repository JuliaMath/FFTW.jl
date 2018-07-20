# This file was formerly a part of Julia. License is MIT: https://julialang.org/license
using FFTW
using FFTW: fftw_vendor
using AbstractFFTs: Plan, plan_inv
using Test
using LinearAlgebra
using Compat

# Base Julia issue #19892
# (test this first to make sure it happens before set_num_threads)
let a = randn(10^5,1), p1 = plan_rfft(a, flags=FFTW.ESTIMATE)
    FFTW.set_num_threads(2)
    p2 = plan_rfft(a, flags=FFTW.ESTIMATE)
    @test p1*a ≈ p2*a
    # make sure threads are actually being used for p2
    # (tests #21163).
    if FFTW.has_sprint_plan
        @test !occursin("dft-thr", string(p1))
        @test occursin("dft-thr", string(p2))
    end
end

# fft
a = rand(8) + im*rand(8)
@test norm(ifft(fft(a)) - a) < 1e-8
@test norm(ifft(fft(a,1),1) - a) < 1e-8
@test norm(ifft(fft(a,[1]),[1]) - a) < 1e-8
@test norm(ifft(fft(a,(1,)),(1,)) - a) < 1e-8
a = rand(-10:10, 8) + im*rand(-10:10, 8)
@test norm(ifft(fft(a)) - a) < 1e-8

m4 = [16.    2     3    13;
    5    11    10     8;
    9     7     6    12;
    4    14    15     1]

true_fft_m4 = [
    34.            34.            34.            34.;
     7. - 1.0im  -5. + 3.0im  -3. + 5.0im   1. - 7.0im;
    16.           -16.           -16.            16.;
     7. + 1.0im  -5. - 3.0im  -3. - 5.0im   1. + 7.0im ]

true_fftn_m4 = [
 136.        0          0         0 ;
   0.       20          8 + 8im   0 - 12im ;
   0.       32 + 32im   0        32 - 32im ;
   0.        0 + 12im   8 - 8im  20 ]

true_fftd2_m4 = [
   34.   13 + 11im    4   13 - 11im ;
   34.   -5 -  3im   -4   -5 +  3im ;
   34.    3 +  5im   -4    3 -  5im ;
   34.  -11 - 13im    4  -11 + 13im ]

b = rand(17,14)
b[3:6,9:12] = m4
sm4 = view(b,3:6,9:12)

m3d = map(Float32,copy(reshape(1:5*3*2, 5, 3, 2)))
true_fftd3_m3d = Array{Float32}(undef, 5, 3, 2)
true_fftd3_m3d[:,:,1] = 17:2:45
true_fftd3_m3d[:,:,2] .= -15

# use invoke to force usage of CTPlan versions even if FFTW is present
for A in (Array,SubArray)
    for f in (:fft,:ifft,:plan_fft,:plan_ifft)
        f_ = Symbol(f, "_")
        @eval begin
            $f_(x::$A{T,N}) where {T,N} = invoke($f, Tuple{AbstractArray{T,N}}, x)
            $f_(x::$A{T,N},r::R) where {T,N,R} = invoke($f,Tuple{AbstractArray{T,N},R},x,r)
        end
    end
end

for (f,fi,pf,pfi) in ((fft,ifft,plan_fft,plan_ifft),
                      (fft_,ifft_,plan_fft_,plan_ifft_))
    pm4 = pf(m4,1)

    fft_m4 = f(m4,1)
    fftd2_m4 = f(m4,2)
    ifft_fft_m4 = fi(f(m4,1),1)
    fftn_m4 = f(m4)
    ifftn_fftn_m4 = fi(f(m4))

    fft!_m4 = complex(m4); fft!(fft!_m4,1)
    fft!d2_m4 = complex(m4); fft!(fft!d2_m4,2)
    ifft!_fft_m4 = f(m4,1); ifft!(ifft!_fft_m4,1)
    fft!n_m4 = complex(m4); fft!(fft!n_m4)
    ifft!n_fftn_m4 = f(m4); ifft!(ifft!n_fftn_m4)

    pfft_m4 = pf(m4,1)*m4
    pfftd2_m4 = pf(m4,2)*m4
    pifft_fft_m4 = pfi(fft_m4,1)*fft_m4
    pfftn_m4 = pf(m4)*m4
    pifftn_fftn_m4 = pfi(fftn_m4)*fftn_m4

    pfft!_m4 = complex(m4); plan_fft!(pfft!_m4,1)*pfft!_m4
    pfft!d2_m4 = complex(m4); plan_fft!(pfft!d2_m4,2)*pfft!d2_m4
    pifft!_fft_m4 = f(m4,1); plan_ifft!(pifft!_fft_m4,1)*pifft!_fft_m4
    pfft!n_m4 = complex(m4); plan_fft!(pfft!n_m4)*pfft!n_m4
    pifft!n_fftn_m4 = f(m4); plan_ifft!(pifft!n_fftn_m4)*pifft!n_fftn_m4

    sfftn_m4 = f(sm4)
    psfftn_m4 = pf(sm4)*sm4
    sfft!n_b = map(Complex{Float64},b)
    sfft!n_m4 = view(sfft!n_b,3:6,9:12); fft!(sfft!n_m4)
    psfft!n_b = map(Complex{Float64},b)
    psfft!n_m4 = view(psfft!n_b,3:6,9:12); plan_fft!(psfft!n_m4)*psfft!n_m4

    for i = 1:length(m4)
        @test fft_m4[i] ≈ true_fft_m4[i]
        @test fftd2_m4[i] ≈ true_fftd2_m4[i]
        @test ifft_fft_m4[i] ≈ m4[i]
        @test fftn_m4[i] ≈ true_fftn_m4[i]
        @test ifftn_fftn_m4[i] ≈ m4[i]

        @test fft!_m4[i] ≈ true_fft_m4[i]
        @test fft!d2_m4[i] ≈ true_fftd2_m4[i]
        @test ifft!_fft_m4[i] ≈ m4[i]
        @test fft!n_m4[i] ≈ true_fftn_m4[i]
        @test ifft!n_fftn_m4[i] ≈ m4[i]

        @test pfft_m4[i] ≈ true_fft_m4[i]
        @test pfftd2_m4[i] ≈ true_fftd2_m4[i]
        @test pifft_fft_m4[i] ≈ m4[i]
        @test pfftn_m4[i] ≈ true_fftn_m4[i]
        @test pifftn_fftn_m4[i] ≈ m4[i]

        @test pfft!_m4[i] ≈ true_fft_m4[i]
        @test pfft!d2_m4[i] ≈ true_fftd2_m4[i]
        @test pifft!_fft_m4[i] ≈ m4[i]
        @test pfft!n_m4[i] ≈ true_fftn_m4[i]
        @test pifft!n_fftn_m4[i] ≈ m4[i]

        @test sfftn_m4[i] ≈ true_fftn_m4[i]
        @test sfft!n_m4[i] ≈ true_fftn_m4[i]
        @test psfftn_m4[i] ≈ true_fftn_m4[i]
        @test psfft!n_m4[i] ≈ true_fftn_m4[i]
    end

    ifft!(sfft!n_m4)
    plan_ifft!(psfft!n_m4)*psfft!n_m4
    @test norm(sfft!n_m4 - m4) < 1e-8
    @test norm(psfft!n_m4 - m4) < 1e-8

    # The following capabilities are FFTW only.
    # They are not available in MKL, and hence do not test them.
    if fftw_vendor != :mkl
        ifft3_fft3_m3d = fi(f(m3d))

        fftd3_m3d = f(m3d,3)
        ifftd3_fftd3_m3d = fi(fftd3_m3d,3)

        fft!d3_m3d = complex(m3d); fft!(fft!d3_m3d,3)
        ifft!d3_fftd3_m3d = copy(fft!d3_m3d); ifft!(ifft!d3_fftd3_m3d,3)

        pfftd3_m3d = pf(m3d,3)*m3d
        pifftd3_fftd3_m3d = pfi(fftd3_m3d,3)*fftd3_m3d

        pfft!d3_m3d = complex(m3d); plan_fft!(pfft!d3_m3d,3)*pfft!d3_m3d
        pifft!d3_fftd3_m3d = copy(fft!d3_m3d); plan_ifft!(pifft!d3_fftd3_m3d,3)*pifft!d3_fftd3_m3d

        @test isa(fftd3_m3d, Array{Complex{Float32},3})
        @test isa(ifftd3_fftd3_m3d, Array{Complex{Float32},3})
        @test isa(fft!d3_m3d, Array{Complex{Float32},3})
        @test isa(ifft!d3_fftd3_m3d, Array{Complex{Float32},3})
        @test isa(pfftd3_m3d, Array{Complex{Float32},3})
        @test isa(pifftd3_fftd3_m3d, Array{Complex{Float32},3})
        @test isa(pfft!d3_m3d, Array{Complex{Float32},3})
        @test isa(pifft!d3_fftd3_m3d, Array{Complex{Float32},3})

        for i = 1:length(m3d)
            @test fftd3_m3d[i] ≈ true_fftd3_m3d[i]
            @test ifftd3_fftd3_m3d[i] ≈ m3d[i]
            @test ifft3_fft3_m3d[i] ≈ m3d[i]

            @test fft!d3_m3d[i] ≈ true_fftd3_m3d[i]
            @test ifft!d3_fftd3_m3d[i] ≈ m3d[i]

            @test pfftd3_m3d[i] ≈ true_fftd3_m3d[i]
            @test pifftd3_fftd3_m3d[i] ≈ m3d[i]
            @test pfft!d3_m3d[i] ≈ true_fftd3_m3d[i]
            @test pifft!d3_fftd3_m3d[i] ≈ m3d[i]
        end
    end  # if fftw_vendor != :mkl

end

# rfft/rfftn

rfft_m4 = rfft(m4,1)
rfftd2_m4 = rfft(m4,2)
rfftn_m4 = rfft(m4)

prfft_m4 = plan_rfft(m4,1)*m4
prfftd2_m4 = plan_rfft(m4,2)*m4
prfftn_m4 = plan_rfft(m4)*m4

srfftn_m4 = rfft(sm4)
psrfftn_m4 = plan_rfft(sm4)*sm4

for i = 1:3, j = 1:4
    @test rfft_m4[i,j] ≈ true_fft_m4[i,j]
    @test rfftd2_m4[j,i] ≈ true_fftd2_m4[j,i]
    @test rfftn_m4[i,j] ≈ true_fftn_m4[i,j]

    @test prfft_m4[i,j] ≈ true_fft_m4[i,j]
    @test prfftd2_m4[j,i] ≈ true_fftd2_m4[j,i]
    @test prfftn_m4[i,j] ≈ true_fftn_m4[i,j]

    @test srfftn_m4[i,j] ≈ true_fftn_m4[i,j]
    @test psrfftn_m4[i,j] ≈ true_fftn_m4[i,j]
end

irfft_rfft_m4 = irfft(rfft_m4,size(m4,1),1)
irfft_rfftd2_m4 = irfft(rfftd2_m4,size(m4,2),2)
irfftn_rfftn_m4 = irfft(rfftn_m4,size(m4,1))

pirfft_rfft_m4 = plan_irfft(rfft_m4,size(m4,1),1)*rfft_m4
pirfft_rfftd2_m4 = plan_irfft(rfftd2_m4,size(m4,2),2)*rfftd2_m4
pirfftn_rfftn_m4 = plan_irfft(rfftn_m4,size(m4,1))*rfftn_m4

for i = 1:length(m4)
    @test irfft_rfft_m4[i] ≈ m4[i]
    @test irfft_rfftd2_m4[i] ≈ m4[i]
    @test irfftn_rfftn_m4[i] ≈ m4[i]

    @test pirfft_rfft_m4[i] ≈ m4[i]
    @test pirfft_rfftd2_m4[i] ≈ m4[i]
    @test pirfftn_rfftn_m4[i] ≈ m4[i]
end

# rfft/rfftn with preallocated array

prfft_m4_prealloc = zero(prfft_m4)
prfftd2_m4_prealloc = zero(prfftd2_m4)
prfftn_m4_prealloc = zero(prfftn_m4)

pirfft_rfft_m4_prealloc = zero(pirfft_rfft_m4)
pirfft_rfftd2_m4_prealloc = zero(pirfft_rfftd2_m4)
pirfftn_rfftn_m4_prealloc = zero(pirfftn_rfftn_m4)

mul!(prfft_m4_prealloc,plan_rfft(m4,1),m4)
mul!(prfftd2_m4_prealloc,plan_rfft(m4,2),m4)
mul!(prfftn_m4_prealloc,plan_rfft(m4),m4)

for i = 1:3, j = 1:4
    @test prfft_m4_prealloc[i,j] ≈ true_fft_m4[i,j]
    @test prfftd2_m4_prealloc[j,i] ≈ true_fftd2_m4[j,i]
    @test prfftn_m4_prealloc[i,j] ≈ true_fftn_m4[i,j]
end

mul!(pirfft_rfft_m4_prealloc,plan_irfft(rfft_m4,size(m4,1),1),prfft_m4_prealloc)
mul!(pirfft_rfftd2_m4_prealloc,plan_irfft(rfftd2_m4,size(m4,2),2),prfftd2_m4_prealloc)
mul!(pirfftn_rfftn_m4_prealloc,plan_irfft(rfftn_m4,size(m4,1)),prfftn_m4_prealloc)

for i = 1:length(m4)
    @test pirfft_rfft_m4_prealloc[i] ≈ m4[i]
    @test pirfft_rfftd2_m4_prealloc[i] ≈ m4[i]
    @test pirfftn_rfftn_m4_prealloc[i] ≈ m4[i]
end

if fftw_vendor != :mkl
    rfftn_m3d = rfft(m3d)
    rfftd3_m3d = rfft(m3d,3)
    @test size(rfftd3_m3d) == size(true_fftd3_m3d)
    irfft_rfftd3_m3d = irfft(rfftd3_m3d,size(m3d,3),3)
    irfftn_rfftn_m3d = irfft(rfftn_m3d,size(m3d,1))
    for i = 1:length(m3d)
        @test rfftd3_m3d[i] ≈ true_fftd3_m3d[i]
        @test irfft_rfftd3_m3d[i] ≈ m3d[i]
        @test irfftn_rfftn_m3d[i] ≈ m3d[i]
    end

    fftn_m3d = fft(m3d)
    @test size(fftn_m3d) == (5,3,2)
    rfftn_m3d = rfft(m3d)
    @test size(rfftn_m3d) == (3,3,2)
    for i = 1:3, j = 1:3, k = 1:2
        @test rfftn_m3d[i,j,k] ≈ fftn_m3d[i,j,k]
    end
end # !mkl

# FFT self-test algorithm (for unscaled 1d forward FFTs):
#   Funda Ergün, "Testing multivariate linear functions: Overcoming
#   the generator bottleneck," Proc. 27th ACM Symposium on the Theory
#   of Computing, pp. 407-416 (1995).
# Check linearity, impulse-response, and time-shift properties.
function fft_test(p::Plan{T}, ntrials=4,
                  tol=1e5 * eps(real(T))) where T<:Complex
    ndims(p) == 1 || throw(ArgumentError("not a 1d FFT"))
    n = length(p)
    twopi_i = (-2 * convert(real(T), π)/n * (0:n-1)) * im
    for trial = 1:ntrials
        # linearity:
        x = rand(T, n)
        y = rand(T, n)
        α = rand(T)
        β = rand(T)
        X = p * (α*x + β*y)
        err = norm(α * (p*x) + β * (p*y) - X, Inf) / norm(X, Inf)
        err <= tol || error("linearity error $err in $p")

        # impulse-response:
        z = zeros(T, n)
        i = rand(0:n-1)
        z[i+1] = 1
        X = exp.(twopi_i*i)
        err = norm(p*z - X, Inf) / norm(X, Inf)
        err <= tol || error("impulse-response error $err in $p")

        # time-shift:
        if n > 1
            s = rand(1:n-1)
            X = (p*x).*exp.(twopi_i*s)
            err = norm(p*circshift(x,s) - X, Inf) / norm(X, Inf)
            err <= tol || error("time-shift error $err in $p")
        end
    end
end

for T in (Complex{Float32}, Complex{Float64})
    for n in [1:100; 121; 143; 1000; 1024; 1031; 2000; 2048]
        x = zeros(T, n)
        fft_test(plan_fft(x))
        fft_test(plan_fft_(x))
    end
end

# test inversion, scaling, and pre-allocated variants
for T in (Complex{Float32}, Complex{Float64})
    for x in (T[1:100;], copy(reshape(T[1:200;], 20,10)))
        y = similar(x)
        for planner in (plan_fft, plan_fft_, plan_ifft, plan_ifft_)
            p = planner(x)
            pi = inv(p)
            p3 = 3*p
            p3i = inv(p3)
            @test eltype(p) == eltype(pi) == eltype(p3) == eltype(p3i) == T
            @test norm(x - p3i * (p * 3x)) < eps(real(T)) * 10000
            @test norm(3x - pi * (p3 * x)) < eps(real(T)) * 10000
            mul!(y, p, x)
            @test y == p * x
            ldiv!(y, p, x)
            @test y == p \ x
        end
    end
end

let
    plan32 = plan_fft([1.0:2048.0;])
    plan64 = plan_fft([1f0:2048f0;])
    FFTW.flops(plan32)
    FFTW.flops(plan64)
end

# Base Julia issue #9772
for x in (randn(10),randn(10,12))
    z = complex(x)
    y = rfft(x)
    @inferred rfft(x)

    # See Julia issue #23063
    if ndims(x) == 2
        @test_broken @inferred brfft(x,18)
    end

    @inferred brfft(y,10)
    for f in (plan_bfft!, plan_fft!, plan_ifft!,
              plan_bfft, plan_fft, plan_ifft,
              fft, bfft, fft_, ifft)
        p = @inferred f(z)
        if isa(p, Plan)
            @inferred plan_inv(p)
        end
    end
    for f in (plan_bfft, plan_fft, plan_ifft,
              plan_rfft, fft, bfft, fft_, ifft)
        # More of #23063 (why does plan_rfft work and the others don't)?
        if ndims(x) == 2 && f != plan_rfft
            @test_broken @inferred f(x)
            @test_broken @inferred plan_inv(f(x))
            continue
        end
        p = @inferred f(x)
        if isa(p, Plan)
            @inferred plan_inv(p)
        end
    end
    # note: inference doesn't work for plan_fft_ since the
    #       algorithm steps are included in the CTPlan type
end

# Base Julia issue #17896
a = rand(5)
@test  fft(a) ==  fft(view(a,:)) ==  fft(view(a, 1:5)) ==  fft(view(a, [1:5;]))
@test rfft(a) == rfft(view(a,:)) == rfft(view(a, 1:5)) == rfft(view(a, [1:5;]))
a16 = convert(Vector{Float16}, a)
@test  fft(a16) ==  fft(view(a16,:)) ==  fft(view(a16, 1:5)) ==  fft(view(a16, [1:5;]))
@test rfft(a16) == rfft(view(a16,:)) == rfft(view(a16, 1:5)) == rfft(view(a16, [1:5;]))

# Discrete cosine transform (DCT) tests

if fftw_vendor != :mkl
    a = rand(8,11) + im*rand(8,11)
    @test norm(idct(dct(a)) - a) < 1e-8

    X = reshape([1,2,7,2,1,5,9,-1,3,4,6,9],3,4)
    Y = rand(17,14)
    Y[3:5,9:12] = X
    sX = view(Y,3:5,9:12)

    true_Xdct = [  13.856406460551018  -3.863239728836245   2.886751345948129  -0.274551994240164; -2.828427124746190  -2.184015211898548  -4.949747468305834   3.966116180118245; 4.898979485566356  -0.194137576915510  -2.857738033247041   2.731723009609389 ]

    true_Xdct_1 = [    5.773502691896258   4.618802153517007   6.350852961085884  10.969655114602890; -4.242640687119286  -2.121320343559643   4.242640687119286  -3.535533905932738; 1.632993161855452   2.041241452319315   5.715476066494083   0.408248290463863 ]

    true_Xdct_2 = [    8.  -3.854030797826254  -3.0  3.761176226848022;
        4.0  -2.071929829606556   4.0  -2.388955165168770; 12.  -0.765366864730179   4.0  -1.847759065022573 ]

    Xdct = dct(X)
    Xdct! = float(X); dct!(Xdct!)
    Xdct_1 = dct(X,1)
    Xdct!_1 = float(X); dct!(Xdct!_1,1)
    Xdct_2 = dct(X,2)
    Xdct!_2 = float(X); dct!(Xdct!_2,2)

    Xidct = idct(true_Xdct)
    Xidct! = copy(true_Xdct); idct!(Xidct!)
    Xidct_1 = idct(true_Xdct_1,1)
    Xidct!_1 = copy(true_Xdct_1); idct!(Xidct!_1,1)
    Xidct_2 = idct(true_Xdct_2,2)
    Xidct!_2 = copy(true_Xdct_2); idct!(Xidct!_2,2)

    pXdct = plan_dct(X)*(X)
    pXdct! = float(X); plan_dct!(pXdct!)*(pXdct!)
    pXdct_1 = plan_dct(X,1)*(X)
    pXdct!_1 = float(X); plan_dct!(pXdct!_1,1)*(pXdct!_1)
    pXdct_2 = plan_dct(X,2)*(X)
    pXdct!_2 = float(X); plan_dct!(pXdct!_2,2)*(pXdct!_2)

    pXidct = plan_idct(true_Xdct)*(true_Xdct)
    pXidct! = copy(true_Xdct); plan_idct!(pXidct!)*(pXidct!)
    pXidct_1 = plan_idct(true_Xdct_1,1)*(true_Xdct_1)
    pXidct!_1 = copy(true_Xdct_1); plan_idct!(pXidct!_1,1)*(pXidct!_1)
    pXidct_2 = plan_idct(true_Xdct_2,2)*(true_Xdct_2)
    pXidct!_2 = copy(true_Xdct_2); plan_idct!(pXidct!_2,2)*(pXidct!_2)

    sXdct = dct(sX)
    psXdct = plan_dct(sX)*(sX)
    sYdct! = copy(Y); sXdct! = view(sYdct!,3:5,9:12); dct!(sXdct!)
    psYdct! = copy(Y); psXdct! = view(psYdct!,3:5,9:12); plan_dct!(psXdct!)*(psXdct!)

    for i = 1:length(X)
        @test Xdct[i] ≈ true_Xdct[i]
        @test Xdct![i] ≈ true_Xdct[i]
        @test Xdct_1[i] ≈ true_Xdct_1[i]
        @test Xdct!_1[i] ≈ true_Xdct_1[i]
        @test Xdct_2[i] ≈ true_Xdct_2[i]
        @test Xdct!_2[i] ≈ true_Xdct_2[i]

        @test pXdct[i] ≈ true_Xdct[i]
        @test pXdct![i] ≈ true_Xdct[i]
        @test pXdct_1[i] ≈ true_Xdct_1[i]
        @test pXdct!_1[i] ≈ true_Xdct_1[i]
        @test pXdct_2[i] ≈ true_Xdct_2[i]
        @test pXdct!_2[i] ≈ true_Xdct_2[i]

        @test Xidct[i] ≈ X[i]
        @test Xidct![i] ≈ X[i]
        @test Xidct_1[i] ≈ X[i]
        @test Xidct!_1[i] ≈ X[i]
        @test Xidct_2[i] ≈ X[i]
        @test Xidct!_2[i] ≈ X[i]

        @test pXidct[i] ≈ X[i]
        @test pXidct![i] ≈ X[i]
        @test pXidct_1[i] ≈ X[i]
        @test pXidct!_1[i] ≈ X[i]
        @test pXidct_2[i] ≈ X[i]
        @test pXidct!_2[i] ≈ X[i]

        @test sXdct[i] ≈ true_Xdct[i]
        @test psXdct[i] ≈ true_Xdct[i]
        @test sXdct![i] ≈ true_Xdct[i]
        @test psXdct![i] ≈ true_Xdct[i]
    end
end # fftw_vendor != :mkl

# test UNALIGNED flag
let A = rand(Float32, 35), Ac = rand(Complex{Float32}, 35)
    local Y = Array{Complex{Float32}}(undef, 20)
    Yc = Array{Complex{Float32}}(undef, 35)
    planr = plan_rfft(Array{Float32}(undef, 32), flags=FFTW.UNALIGNED)
    planc = plan_fft(Array{Complex{Float32}}(undef, 32), flags=FFTW.UNALIGNED)
    for ioff in 0:3
        ii = 1+ioff:32+ioff
        @test planr * view(A, ii) ≈ planr * A[ii] ≈ rfft(view(A, ii)) ≈ rfft(A[ii])
        @test planc * view(Ac, ii) ≈ planc * Ac[ii] ≈ fft(view(Ac, ii)) ≈ fft(Ac[ii])
        for ooff in 0:3
            io = 1+ooff:17+ooff
            FFTW.mul!(view(Y, io), planr, view(A, ii))
            @test Y[io] ≈ rfft(A[ii])
            FFTW.mul!(view(Y, io), planr, A[ii])
            @test Y[io] ≈ rfft(A[ii])
            io = 1+ooff:32+ooff
            FFTW.mul!(view(Yc, io), planc, view(Ac, ii))
            @test Yc[io] ≈ fft(Ac[ii])
            FFTW.mul!(view(Yc, io), planc, Ac[ii])
            @test Yc[io] ≈ fft(Ac[ii])
        end
    end
    @test_throws ArgumentError plan_rfft(Array{Float32}(undef, 32)) * view(A, 2:33)
    @test_throws ArgumentError plan_fft(Array{Complex{Float32}}(undef, 32)) * view(Ac, 2:33)
end
