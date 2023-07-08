module FFTWChainRulesCoreExt

using FFTW
using FFTW: r2r
using ChainRulesCore

# DCT

function ChainRulesCore.frule(Δ, ::typeof(dct), x::AbstractArray, region = 1:ndims(x))
    Δx = Δ[2]
    y = dct(x, region)
    Δy = dct(Δx, region)
    return y, Δy
end

function ChainRulesCore.rrule(::typeof(dct), x::AbstractArray, region...)
    y = dct(x, region...)
    project_x = ProjectTo(x)

    function dct_pullback(ȳ)
        f̄ = NoTangent()
        x̄ = project_x(idct(unthunk(ȳ), region...))
        r̄ = NoTangent()

        if isempty(region)
            return f̄, x̄
        else
            return f̄, x̄, r̄
        end
    end

    return y, dct_pullback
end

# IDCT

function ChainRulesCore.frule(Δ, ::typeof(idct), x::AbstractArray, region = 1:ndims(x))
    Δx = Δ[2]
    y = idct(x, region)
    Δy = idct(Δx, region)
    return y, Δy
end

function ChainRulesCore.rrule(::typeof(idct), x::AbstractArray, region...)
    y = idct(x, region...)
    project_x = ProjectTo(x)

    function idct_pullback(ȳ)
        f̄ = NoTangent()
        x̄ = project_x(dct(unthunk(ȳ), region...))
        r̄ = NoTangent()

        if isempty(region)
            return f̄, x̄
        else
            return f̄, x̄, r̄
        end
    end

    return y, idct_pullback
end

# R2R

function ChainRulesCore.frule(Δ, ::typeof(r2r), x::AbstractArray, kind, region = 1:ndims(x))
    Δx = Δ[2]
    y = r2r(x, kind, region)
    Δy = r2r(Δx, kind, region)
    return y, Δy
end

end # module
