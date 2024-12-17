module FFTWForwardDiffExt
using FFTW
using ForwardDiff
import FFTW: plan_r2r, r2r
import FFTW.AbstractFFTs: dualplan, dual2array
import ForwardDiff: Dual

plan_r2r(x::AbstractArray{D}, FLAG, dims=1:ndims(x)) where D<:Dual = dualplan(D, plan_r2r(dual2array(x), FLAG, 1 .+ dims))
plan_r2r(x::AbstractArray{<:Complex{D}}, FLAG, dims=1:ndims(x)) where D<:Dual = dualplan(D, plan_r2r(dual2array(x), FLAG, 1 .+ dims))


r2r(x::AbstractArray{<:Dual}, kinds, region...) = plan_r2r(x, kinds, region...) * x
r2r(x::AbstractArray{<:Complex{<:Dual}}, kinds, region...) = plan_r2r(x, kinds, region...) * x

end #module