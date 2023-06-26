module FFTWChainRulesCoreExt

using FFTW
using FFTW: r2r
using ChainRulesCore

# DCT

function ChainRulesCore.frule(Δ, ::typeof(dct), x::AbstractArray, args...)
    Δx = Δ[2]
    y = dct(x, args...)
    Δy = dct(Δx, args...)
    return y, Δy
end

function ChainRulesCore.rrule(::typeof(dct), x::AbstractArray, args...)
    y = dct(x, args...)
    project_x = ProjectTo(x)

    function dct_pullback(ȳ)
        f̄ = NoTangent()
        x̄ = project_x(idct(unthunk(ȳ), args...))
        ā = NoTangent()

        if isempty(args)
            return f̄, x̄
        else
            return f̄, x̄, ā
        end
    end

    return y, dct_pullback
end

# IDCT

function ChainRulesCore.frule(Δ, ::typeof(idct), x::AbstractArray, args...)
    Δx = Δ[2]
    y = idct(x, args...)
    Δy = idct(Δx, args...)
    return y, Δy
end

function ChainRulesCore.rrule(::typeof(idct), x::AbstractArray, args...)
    y = idct(x, args...)
    project_x = ProjectTo(x)

    function idct_pullback(ȳ)
        f̄ = NoTangent()
        x̄ = project_x(dct(unthunk(ȳ), args...))
        ā = NoTangent()

        if isempty(args)
            return f̄, x̄
        else
            return f̄, x̄, ā
        end
    end

    return y, idct_pullback
end

# R2R

function ChainRulesCore.frule(Δ, ::typeof(r2r), x::AbstractArray, args...)
    Δx = Δ[2]
    y = r2r(x, args...)
    Δy = r2r(Δx, args...)
    return y, Δy
end

function ChainRulesCore.rrule(::typeof(r2r), x::AbstractArray, kinds, args...)
    y = r2r(x, kinds, args...)
    kinvs = Tuple(FFTW.inv_kind[k] for k in kinds)
    project_x = ProjectTo(x)

    function r2r_pullback(ȳ)
        f̄ = NoTangent()
        x̄ = project_x(r2r(unthunk(ȳ), kinvs, args...))
        k̄ = NoTangent()
        ā = NoTangent()

        if isempty(args)
            return f̄, x̄, k̄
        else
            return f̄, x̄, k̄, ā
        end
    end

    return y, r2r_pullback
end

end # module
