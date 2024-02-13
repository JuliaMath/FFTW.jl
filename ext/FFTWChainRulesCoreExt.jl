module FFTWChainRulesCoreExt

using FFTW
using FFTW: r2r
using ChainRulesCore

# DCT/IDCT

for (fwd, bwd) in (
    (dct, idct),
    (idct, dct),
)
    function ChainRulesCore.frule(Δ, ::typeof(fwd), x::AbstractArray, region = 1:ndims(x))
        Δx = Δ[2]
        y = fwd(x, region)
        Δy = fwd(Δx, region)
        return y, Δy
    end

    function ChainRulesCore.rrule(::typeof(fwd), x::AbstractArray)
        project_x = ProjectTo(x)
        dct_pb(Δ) = NoTangent(), project_x(bwd(unthunk(Δ)))
        return fwd(x), dct_pb
    end

    function ChainRulesCore.rrule(::typeof(fwd), x::AbstractArray, region)
        project_x = ProjectTo(x)
        dct_pb(Δ) = NoTangent(), project_x(bwd(unthunk(Δ), region)), NoTangent()
        return fwd(x, region), dct_pb
    end
end

# R2R

function ChainRulesCore.frule(Δ, ::typeof(r2r), x::AbstractArray, kind, region = 1:ndims(x))
    Δx = Δ[2]
    y = r2r(x, kind, region)
    Δy = r2r(Δx, kind, region)
    return y, Δy
end

end # module
