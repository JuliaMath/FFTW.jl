const IODIM = Tuple{Int,Tuple{Int,Int}}
might_reshape!(iodims::Vector{IODIM}) = begin
    len = 1
    @inbounds for j in 2:length(iodims)
        szʲ, stsʲ = iodims[j]
        i = len
        while i > 0
            szⁱ, stsⁱ = iodims[i]
            if stsⁱ .* szⁱ == stsʲ # able to reshape
                iodims[i] = szⁱ * szʲ, stsⁱ
                break
            end
            i -= 1
        end
        i > 0 && continue
        iodims[len+=1] = szʲ, stsʲ # fail to reshape
    end
    resize!(iodims, len)
end

howmany_loop(iodims::Vector{IODIM}) = begin
    length(iodims) <= 1 && return iodims, IODIM[]
    iodims = sort!(iodims, by = x -> x[2][1]) |> might_reshape!
    sort!(iodims, by = x -> (x[2], -x[1]))
    howmany, loop = [popfirst!(iodims)], iodims
    @inbounds for i in length(loop):-1:2
        (szʲ, stsʲ), (szⁱ, stsⁱ) = loop[i-1], loop[i]
        loop[i] = szⁱ, stsⁱ .- szʲ .* stsʲ
    end
    howmany, loop
end

dims_howmany_loop(X::StridedArray, Y::StridedArray, sz, region) = begin
    reg = unique(region)
    length(reg) != length(region) &&
        throw(ArgumentError("each dimension can be transformed at most once"))
    iodims = [(x,(y,z)) for (x,y,z) in zip(sz, strides(X), strides(Y))]
    oreg = [1:ndims(X)...]
    oreg[reg] .= 0
    oreg .*= sz .> 1
    oreg = findall(>(0), oreg)
    dims = iodims[reg]
    howmany, loop = howmany_loop(iodims[oreg])
    dims, howmany, loop
end

# Although MKL's plan has no alignment limitation,
# but the check is preseved here to throw the same error
mutable struct cLoopPlan{T,K,inplace,N,G} <: FFTWPlan{T,K,inplace}
    plan::PlanPtr
    sz::NTuple{N,Int} # size of array on which plan operates (Int tuple)
    osz::NTuple{N,Int} # size of output array (Int tuple)
    istride::NTuple{N,Int} # strides of input
    ostride::NTuple{N,Int} # strides of output
    ialign::Int32 # alignment mod 16 of input
    oalign::Int32 # alignment mod 16 of input
    flags::UInt32 # planner flags
    region::G # region (iterable) of dims that are transormed
    loop::Vector{IODIM} # keep the loop size if needed
    pinv::ScaledPlan
    function cLoopPlan{T,K,inplace}(plan::PlanPtr, X, Y, flags, region, loop) where {T,K,inplace}
        N, G = ndims(X), typeof(region)
        p = new{T,K,inplace,N,G}(plan, size(X), size(Y), strides(X), strides(Y),
                                    alignment_of(X), alignment_of(Y), flags,
                                    region, loop)
        finalizer(maybe_destroy_plan, p)
        p
    end
end

for (Tr,Tc,fftw,lib) in ((:Float64,:(Complex{Float64}),"fftw",:libfftw3),
                         (:Float32,:(Complex{Float32}),"fftwf",:libfftw3f))
    @eval @exclusive function cLoopPlan{$Tc,K,inplace}(X::StridedArray{$Tc,N},
                                        Y::StridedArray{$Tc,N}, region, flags, timelimit) where {K,inplace,N}
        unsafe_set_timelimit($Tr, timelimit)
        dims, howmany, loop = dims_howmany_loop(X, Y, size(X), region)
        plan = ccall(($(string(fftw,"_plan_guru64_dft")),$lib),
                     PlanPtr,
                     (Int32, Ptr{Int}, Int32, Ptr{Int},
                      Ptr{$Tc}, Ptr{$Tc}, Int32, UInt32),
                     length(dims), dims, length(howmany), howmany,
                     X, Y, K, UNALIGNED) ## flags is useless
        unsafe_set_timelimit($Tr, NO_TIMELIMIT)
        if plan == C_NULL
            error("FFTW could not create plan") # shouldn't normally happen
        end
        return cLoopPlan{$Tc,K,inplace}(plan, X, Y, flags, region, loop)
    end
end

show(io::IO, p::cLoopPlan{T,K,inplace}) where {T,K,inplace} = begin
    print(io, inplace ? "FFTW in-place " : "FFTW ",
          K < 0 ? "forward" : "backward", " plan for ")
    showfftdims(io, p.sz, p.istride, T)
    has_sprint_plan && print(io, "\n", sprint_plan(p))
end

@generated function unsafe_nd_execute!(p::cLoopPlan{T}, X::Ptr{T}, Y::Ptr{T}, ::Val{N}) where {T,N}
    lib = T <: fftwSingle ? :libfftw3f : :libfftw3
    fun = T <: fftwSingle ? :fftwf_execute_dft : :fftw_execute_dft
    ker = :(@ccall $lib.$fun(p::PlanPtr, X::Ptr{T}, Y::Ptr{T})::Cvoid)
    for i in Base.OneTo(N)
        itr, sz, sts = Symbol(:i,i), Symbol(:sz,i), Symbol(:sts,i)
        ker = :(for $itr in Base.OneTo($sz)
                $ker
                X, Y = (X, Y) .+ $sts
            end
        )
    end
    ex = Expr(:block)
    for i in Base.OneTo(N)
        sz, sts = Symbol(:sz,i), Symbol(:sts,i)
        push!(ex.args, :(@inbounds $sz, $sts = p.loop[$i]))
        push!(ex.args, :($sts = $sts .* $(sizeof(T))))
    end
    push!(ex.args, ker)
    return ex
end

function unsafe_execute!(p::cLoopPlan{T}, x::StridedArray{T},
                         y::StridedArray{T}) where {T <: fftwComplex}
    X, Y = pointer(x), pointer(y)
    loop = p.loop
    length(loop) == 0 && return unsafe_nd_execute!(p, X, Y, Val(0))
    length(loop) == 1 && return unsafe_nd_execute!(p, X, Y, Val(1))
    length(loop) == 2 && return unsafe_nd_execute!(p, X, Y, Val(2))
    length(loop) == 3 && return unsafe_nd_execute!(p, X, Y, Val(3))
    length(loop) == 4 && return unsafe_nd_execute!(p, X, Y, Val(4))
    length(loop) == 5 && return unsafe_nd_execute!(p, X, Y, Val(5))
    unsafe_nd_execute!(p, X, Y, Val(length(loop)))
end

function mul!(y::StridedArray{T,N}, p::cLoopPlan{T}, x::StridedArray{T,N}) where {T,N}
    assert_applicable(p, x, y)
    unsafe_execute!(p, x, y)
    y
end

function *(p::cLoopPlan{T,K,inplace}, x::StridedArray{T,N}) where {T,K,N,inplace}
    assert_applicable(p, x)
    y = inplace ? x : Array{T,N}(undef, p.osz)
    unsafe_execute!(p, x, y)
    y
end

MightNeedLoop{T} = Union{StridedArray{T,3}, StridedArray{T,4}, StridedArray{T,5}, StridedArray{T,6}, StridedArray{T,7}}
for inplace in (false,true)
    for (f,direction) in ((:fft,FORWARD), (:bfft,BACKWARD))
        plan_f = inplace ? Symbol("plan_",f,"!") : Symbol("plan_",f)
        @eval $plan_f(X::MightNeedLoop{<:fftwComplex}, region;
                    flags::Integer=ESTIMATE,
                    timelimit::Real=NO_TIMELIMIT) = begin
            T, N = eltype(X), ndims(X)
            Y = $inplace ? X : FakeArray{T}(size(X))
            cLoopPlan{T,$direction,$inplace}(X, Y, region, flags, timelimit)
        end
        idirection = -direction
        @eval function plan_inv(p::cLoopPlan{T,$direction,$inplace}) where {T<:fftwComplex}
            X = Array{T}(undef, p.sz)
            Y = $inplace ? X : FakeArray{T}(size(X))
            ScaledPlan(cLoopPlan{T,$idirection,$inplace}(X, Y, p.region,
                                                          p.flags, NO_TIMELIMIT),
                       normalization(X, p.region))
        end
    end
end
