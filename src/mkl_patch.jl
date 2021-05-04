dims_howmany_for_loopplan(X::StridedArray, Y::StridedArray, plandims::Tuple{Vararg{Int}}) = begin
    P = unique(plandims)
    length(P) != length(plandims) &&
        throw(ArgumentError("each dimension can be transformed at most once"))
    reg, oreg  = P[1:end-1], P[end]
    sz, ist, ost = collect.((size(X), strides(X), strides(Y)))
    dims = Matrix(transpose([sz[reg] ist[reg] ost[reg]]))
    howmany = Matrix(transpose([sz[oreg] ist[oreg] ost[oreg]]))
    dims, howmany
end
# Although MKL's plan has no alignment limitation,
# but the check is preseved here to keep similar behavior(throw the same error)
mutable struct cLoopPlan{T,K,inplace,N,G,NW} <: FFTWPlan{T,K,inplace}
    plan::PlanPtr
    sz::NTuple{N,Int} # size of array on which plan operates (Int tuple)
    osz::NTuple{N,Int} # size of output array (Int tuple)
    istride::NTuple{N,Int} # strides of input
    ostride::NTuple{N,Int} # strides of output
    ialign::Int32 # alignment mod 16 of input
    oalign::Int32 # alignment mod 16 of input
    flags::UInt32 # planner flags
    region::G # region (iterable) of dims that are transormed
    loopdims::NTuple{NW,Int}
    pinv::ScaledPlan
    function cLoopPlan{T,K,inplace}(plan::PlanPtr, X, Y, region, flags, loopdims) where {T,K,inplace}
        N, NW, G = ndims(X), length(loopdims), typeof(region)
        p = new{T,K,inplace,N,G,NW}(plan, size(X), size(Y), strides(X), strides(Y),
                                    alignment_of(X), alignment_of(Y), flags, region, loopdims)
        finalizer(maybe_destroy_plan, p)
        p
    end
end

for (Tr,Tc,fftw,lib) in ((:Float64,:(Complex{Float64}),"fftw",:libfftw3),
                         (:Float32,:(Complex{Float32}),"fftwf",:libfftw3f))
    @eval @exclusive function cLoopPlan{$Tc,K,inplace}(X::StridedArray{$Tc,N},
                                        Y::StridedArray{$Tc,N}, region, flags, timelimit) where {K,inplace,N}
        unsafe_set_timelimit($Tr, timelimit)
        plandims, loopdims = split_dim(X, region)
        dims, howmany = dims_howmany_for_loopplan(X, Y, plandims)
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
        return cLoopPlan{$Tc,K,inplace}(plan, X, Y, region, flags, loopdims)
    end
end

function split_dim(X, region)
    dims = filter(i -> !in(i, region), 1:ndims(X))
    # Determine which dimension to be planned by MKL's Dfti
    # I believe the first dimension shoule be selected anyhow for columnmajor layout.
    # Other dimensions might be selected based on size for better thread performance.
    id = !in(1,region) || size(X,dims[1]) > 10Threads.nthreads() ?
         1 : argmax(map(i -> size(X,i), dims))
    loopdims = dims[[1:id-1;id+1:end]] |> NTuple{ndims(X) - length(region) - 1, Int}
    plandims = (region..., dims[id])
    plandims, loopdims
end
show(io::IO, p::cLoopPlan{T,K,inplace}) where {T,K,inplace} = begin
    print(io, inplace ? "FFTW in-place " : "FFTW ",
          K < 0 ? "forward" : "backward", " plan for ")
    showfftdims(io, p.sz, p.istride, T)
    has_sprint_plan && print(io, "\n", sprint_plan(p))
end

unsafe_loop_execute!(plan::cLoopPlan{T}, X::Ptr{T}, Y::Ptr{T}) where {T<:fftwComplex} = begin
    if T <:fftwSingle
        @ccall libfftw3f.fftwf_execute_dft(plan::PlanPtr, X::Ptr{T}, Y::Ptr{T})::Cvoid
    else
        @ccall libfftw3.fftw_execute_dft(plan::PlanPtr, X::Ptr{T}, Y::Ptr{T})::Cvoid
    end
end

function unsafe_execute!(p::cLoopPlan{T}, x::StridedArray{T},
                         y::StridedArray{T}) where {T <: fftwComplex}
    loopistride = map(i -> p.istride[i], p.loopdims)
    loopostride = map(i -> p.ostride[i], p.loopdims)
    ax = map(i -> axes(x,i) .- first(axes(x,i)), p.loopdims)
    Elsz, pointerˣ, pointerʸ = sizeof(T), pointer(x), pointer(y)
    for ind in Iterators.product(ax...)
        pointerˣ′ = pointerˣ + Elsz * sum(loopistride .* ind)
        pointerʸ′ = pointerʸ + Elsz * sum(loopostride .* ind)
        unsafe_loop_execute!(p, pointerˣ′, pointerʸ′)
    end
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
        @eval $plan_f(X::MightNeedLoop{<:fftwComplex}, region; kws...) = $plan_f(X, Tuple(region); kws...)
        # length information is needed for type stability, so I use ntuple here.
        @eval $plan_f(X::MightNeedLoop{<:fftwComplex}; kws...) = $plan_f(X, ntuple(identity, ndims(X)); kws...)
        @eval $plan_f(X::MightNeedLoop{<:fftwComplex}, region::Union{Integer,Tuple};
                    flags::Integer=ESTIMATE,
                    timelimit::Real=NO_TIMELIMIT) = begin
            T, N = eltype(X), ndims(X)
            Y = $inplace ? X : FakeArray{T}(size(X))
            N <= length(region) + 1 && return cFFTWPlan{T,$direction,$inplace,N}(X, Y, region, flags, timelimit)
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
