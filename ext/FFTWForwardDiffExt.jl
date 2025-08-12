module FFTWForwardDiffExt
using FFTW
using ForwardDiff
import FFTW: plan_r2r, plan_r2r!, plan_dct, plan_dct!, plan_idct, plan_idct!, r2r, r2r!, dct, dct!, idct, idct!, fftwReal, REDFT10, REDFT01
import FFTW.AbstractFFTs: dualplan, dual2array
import ForwardDiff: Dual


for plan in (:plan_r2r, :plan_r2r!)
    @eval begin
        $plan(x::AbstractArray{D}, FLAG, dims=1:ndims(x)) where D<:Dual = dualplan(D, $plan(dual2array(x), FLAG, 1 .+ dims))
        $plan(x::AbstractArray{<:Complex{D}}, FLAG, dims=1:ndims(x)) where D<:Dual = dualplan(D, $plan(dual2array(x), FLAG, 1 .+ dims))
    end
end

for f in (:r2r, :r2r!)
    pf = Symbol("plan_", f)
    @eval begin
        $f(x::AbstractArray{<:Dual}, kinds, region...) = $pf(x, kinds, region...) * x
        $f(x::AbstractArray{<:Complex{<:Dual}}, kinds, region...) = $pf(x, kinds, region...) * x
    end
end


for f in (:dct, :dct!, :idct, :idct!)
    pf = Symbol("plan_", f)
    @eval begin
        $f(x::AbstractArray{<:Dual}) = $pf(x) * x
        $f(x::AbstractArray{<:Dual}, region) = $pf(x, region) * x
    end
end

for plan in (:plan_dct, :plan_dct!, :plan_idct, :plan_idct!)
    @eval begin
        $plan(x::AbstractArray{D}, dims=1:ndims(x); kwds...) where D<:Dual = dualplan(D, $plan(dual2array(x), 1 .+ dims; kwds...))
        $plan(x::AbstractArray{<:Complex{D}}, dims=1:ndims(x); kwds...) where D<:Dual = dualplan(D, $plan(dual2array(x), 1 .+ dims; kwds...))
    end
end

end #module