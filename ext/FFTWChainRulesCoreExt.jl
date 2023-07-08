module FFTWChainRulesCoreExt

using FFTW
using FFTW: r2r
using ChainRulesCore

# DCT/IDCT

function ChainRulesCore.frule(Δ, ::typeof(dct), x::AbstractArray, region = 1:ndims(x))
    Δx = Δ[2]
    y = dct(x, region)
    Δy = dct(Δx, region)
    return y, Δy
end

function ChainRulesCore.rrule(::typeof(dct), x::AbstractArray)
    project_x = ProjectTo(x)
    region = 1:ndims(x)
    dct_pb(Δ) = NoTangent(), project_x(idct(unthunk(Δ), region))
    return dct(x, region), dct_pb
end

function ChainRulesCore.rrule(::typeof(dct), x::AbstractArray, region)
    project_x = ProjectTo(x)
    dct_pb(Δ) = NoTangent(), project_x(idct(unthunk(Δ), region)), NoTangent()
    return dct(x, region), dct_pb
end

# IDCT

function ChainRulesCore.frule(Δ, ::typeof(idct), x::AbstractArray, region = 1:ndims(x))
    Δx = Δ[2]
    y = idct(x, region)
    Δy = idct(Δx, region)
    return y, Δy
end

function ChainRulesCore.rrule(::typeof(idct), x::AbstractArray)
    project_x = ProjectTo(x)
    region = 1:ndims(x)
    dct_pb(Δ) = NoTangent(), project_x(dct(unthunk(Δ), region))
    return idct(x, region), dct_pb
end

function ChainRulesCore.rrule(::typeof(idct), x::AbstractArray, region)
    project_x = ProjectTo(x)
    dct_pb(Δ) = NoTangent(), project_x(dct(unthunk(Δ), region)), NoTangent()
    return idct(x, region), dct_pb
end

# R2R

function ChainRulesCore.frule(Δ, ::typeof(r2r), x::AbstractArray, kind, region = 1:ndims(x))
    Δx = Δ[2]
    y = r2r(x, kind, region)
    Δy = r2r(Δx, kind, region)
    return y, Δy
end

end # module
