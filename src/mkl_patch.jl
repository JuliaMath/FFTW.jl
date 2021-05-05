might_reshape!(sz::Vector{Int}, ist::Vector{Int}, ost::Vector{Int}) = begin
    ax = eachindex(sz)
    @inbounds for i in ax, j in ax
        # fix for sz = (10,10) ist = (2,20) ost = (1,10)
        if sz[i] > 1 && sz[j] > 1 && (ist[i], ost[i]) .* sz[i] == (ist[j], ost[j])
            sz[i], sz[j] = sz[i] * sz[j], 1
            return might_reshape!(sz, ist, ost)
        end
    end
    p = sz .> 1
    sz[p], ist[p], ost[p]
end

howmany_loopinfo(sz::Vector{Int}, ist::Vector{Int}, ost::Vector{Int}) = begin
    szʳ, istʳ, ostʳ = might_reshape!(sz, ist, ost) # try to reshape to reduce loop dims
    ind = sortperm(tuple.(istʳ, ostʳ, .-szʳ))
    szˢ, istˢ, ostˢ = szʳ[ind], istʳ[ind], ostʳ[ind]
    pick = first(istˢ) == 1 ? 1 : findfirst(>(10Threads.nthreads()), szˢ) #how to improve?
    howmany = [szˢ[pick], istˢ[pick], ostˢ[pick]]
    rest = [1:pick-1; pick+1:length(szˢ)]
    loopinfo = szˢ[rest], istˢ[rest], ostˢ[rest]
    howmany, loopinfo
end

dims_howmany_loopinfo(X::StridedArray, Y::StridedArray, region) = begin
    reg = unique(region)
    length(reg) != length(region) &&
        throw(ArgumentError("each dimension can be transformed at most once"))
    oreg = filter(!in(reg), 1:ndims(X))
    sz, ist, ost = collect.((size(X), strides(X), strides(Y)))
    # remove dimension with size == 1
    sz₁dim = findall(==(1), sz)
    if length(sz₁dim) > 0
        filter!(!in(sz₁dim), reg)
        filter!(!in(sz₁dim), oreg)
    end
    dims = Matrix(transpose([sz[reg] ist[reg] ost[reg]]))
    howmany, loopinfo = if length(oreg) < 2
        Vector([sz[oreg];ist[oreg];ost[oreg]]), (Int[], Int[], Int[])
    else
        howmany_loopinfo(sz[oreg], ist[oreg], ost[oreg])
    end
    dims, howmany, loopinfo
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
    loopsz::Vector{Int} # keep the loop size if needed
    loopistride::Vector{Int} # keep the loop istride if needed
    loopostride::Vector{Int} # keep the loop ostride if needed
    pinv::ScaledPlan
    function cLoopPlan{T,K,inplace}(plan::PlanPtr, X, Y, flags, region, loopinfo) where {T,K,inplace}
        N, G = ndims(X), typeof(region)
        loopsz, loopistride, loopostride = loopinfo
        loopistride .*= sizeof(T)
        loopostride .*= sizeof(T)
        p = new{T,K,inplace,N,G}(plan, size(X), size(Y), strides(X), strides(Y),
                                    alignment_of(X), alignment_of(Y), flags,
                                    region, loopsz, loopistride, loopostride)
        finalizer(maybe_destroy_plan, p)
        p
    end
end

for (Tr,Tc,fftw,lib) in ((:Float64,:(Complex{Float64}),"fftw",:libfftw3),
                         (:Float32,:(Complex{Float32}),"fftwf",:libfftw3f))
    @eval @exclusive function cLoopPlan{$Tc,K,inplace}(X::StridedArray{$Tc,N},
                                        Y::StridedArray{$Tc,N}, region, flags, timelimit) where {K,inplace,N}
        unsafe_set_timelimit($Tr, timelimit)
        dims, howmany, loopinfo = dims_howmany_loopinfo(X, Y, region)
        plan = ccall(($(string(fftw,"_plan_guru64_dft")),$lib),
                     PlanPtr,
                     (Int32, Ptr{Int}, Int32, Ptr{Int},
                      Ptr{$Tc}, Ptr{$Tc}, Int32, UInt32),
                     size(dims,2), dims, size(howmany,2), howmany,
                     X, Y, K, UNALIGNED) ## flags is useless
        unsafe_set_timelimit($Tr, NO_TIMELIMIT)
        if plan == C_NULL
            error("FFTW could not create plan") # shouldn't normally happen
        end
        return cLoopPlan{$Tc,K,inplace}(plan, X, Y, flags, region, loopinfo)
    end
end

show(io::IO, p::cLoopPlan{T,K,inplace}) where {T,K,inplace} = begin
    print(io, inplace ? "FFTW in-place " : "FFTW ",
          K < 0 ? "forward" : "backward", " plan for ")
    showfftdims(io, p.sz, p.istride, T)
    has_sprint_plan && print(io, "\n", sprint_plan(p))
end

unsafe_single_execute!(p::cLoopPlan{T}, X::Ptr{T}, Y::Ptr{T}) where {T<:fftwComplex} = begin
    if T <:fftwSingle
        @ccall libfftw3f.fftwf_execute_dft(p::PlanPtr, X::Ptr{T}, Y::Ptr{T})::Cvoid
    else
        @ccall libfftw3.fftw_execute_dft(p::PlanPtr, X::Ptr{T}, Y::Ptr{T})::Cvoid
    end
end

function unsafe_nd_execute!(p::cLoopPlan{T}, X::Ptr{T}, Y::Ptr{T},
                            sz, ist, ost) where T
    offset(x, y) = length(x) > 1 ? x[1] * y[1] + offset(Base.tail(x), Base.tail(y)) :
                                   x[1] * y[1]
    for ind in Iterators.product(map(sz -> (0:sz-1), sz)...)
        X′ = X + offset(ist, ind)
        Y′ = Y + offset(ost, ind)
        unsafe_single_execute!(p, X′, Y′)
    end
end

function unsafe_execute!(p::cLoopPlan{T}, x::StridedArray{T},
                         y::StridedArray{T}) where {T <: fftwComplex}
    X, Y = pointer(x), pointer(y)
    sz = p.loopsz
    length(sz) == 0 && return unsafe_single_execute!(p, X, Y)
    ist, ost = p.loopistride, p.loopostride
    length(sz) == 1 && return unsafe_nd_execute!(p, X, Y, (sz[1],), (ist[1],), (ost[1],))
    unsafe_nd_execute!(p, X, Y, (sz...,), (ist...,), (ost...,))
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
