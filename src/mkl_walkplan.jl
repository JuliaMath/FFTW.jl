struct PickDim{T,N,P<:StridedArray{T}} <: DenseArray{T,N}
    data::P
    dims::NTuple{N,Int}
    PickDim(data::StridedArray{T}, dims::NTuple{N,Int}) where{T,N} = begin
        new{T,N,typeof(data)}(data, dims)
    end
end
Base.size(A::PickDim) = map(i -> size(A.data)[i], A.dims)
Base.strides(A::PickDim) = map(i -> strides(A.data)[i], A.dims)

mutable struct cWalkPlan{T,K,inplace,N,NW} <: FFTWPlan{T,K,inplace}
    plan::PlanPtr
    sz::NTuple{N,Int}
    osz::NTuple{N,Int}
    istride::NTuple{N,Int}
    ostride::NTuple{N,Int}
    walkdims::NTuple{NW,Int}
    function cWalkPlan{T,K,inplace}(plan::PlanPtr, X, Y, walkdims) where {T,K,inplace}
        p = new{T,K,inplace,ndims(X),length(walkdims)}(plan,
            size(X), size(Y), strides(X), strides(Y), walkdims)
        finalizer(maybe_destroy_plan, p)
        p
    end
end
for (Tr,Tc,fftw,lib) in ((:Float64,:(Complex{Float64}),"fftw",:libfftw3),
                         (:Float32,:(Complex{Float32}),"fftwf",:libfftw3f))
    @eval @exclusive function cWalkPlan{$Tc,direction,inplace}(X::StridedArray{$Tc,N},
                                              Y::StridedArray{$Tc,N},
                                              region, timelimit::Real) where {direction,inplace,N}
        length(region) == length(unique(region)) ||
            throw(ArgumentError("each dimension can be transformed at most once"))
        plandims, walkdims = split_dim(X, region)
        X′, Y′ = PickDim(X, plandims), PickDim(Y, plandims)
        dims, howmany = dims_howmany(X′, Y′, [size(X′)...], 1:length(region))
        plan = ccall(($(string(fftw,"_plan_guru64_dft")),$lib),
                     PlanPtr,
                     (Int32, Ptr{Int}, Int32, Ptr{Int},
                      Ptr{$Tc}, Ptr{$Tc}, Int32, UInt32),
                     size(dims,2), dims, size(howmany,2), howmany,
                     X, Y, direction, 0)
        if plan == C_NULL
            error("FFTW could not create plan") # shouldn't normally happen
        end
        return cWalkPlan{$Tc,K,inplace}(plan, X, Y, walkdims)
    end
end

function split_dim(X, region)
    dims = filter(i -> !in(i, region), 1:ndims(X))
    # Determine which dimension to be planned by MKL's Dfti
    # I believe the first dimension shoule be selected anyhow for columnmajor layout.
    # Other dimensions might be selected based on size for better thread performance.
    id = !in(1,region) || size(X,dims[1]) > 10Threads.nthreads() ?
         1 : argmax(map(i -> size(X,i), dims))
    walkdims = dims[[1:id-1;id+1:end]] |> NTuple{ndims(X) - length(region) - 1, Int}
    plandims = (region...,dims[id])
    plandims, walkdims
end

show(io::IO, p::cWalkPlan{T,K,inplace}) where {T,K,inplace} = begin
    print(io, inplace ? "FFTW in-place " : "FFTW ",
          K < 0 ? "forward" : "backward", " plan for ")
    showfftdims(io, p.sz, p.istride, T)
    has_sprint_plan && print(io, "\n", sprint_plan(p))
end

unsafe_walk_execute!(plan::cWalkPlan{T}, X::Ptr{T}, Y::Ptr{T}) where {T<:fftwComplex} = begin
    if T <:fftwSingle
        @ccall libfftw3f.fftwf_execute_dft(plan::PlanPtr, X::Ptr{T}, Y::Ptr{T})::Cvoid
    else
        @ccall libfftw3.fftw_execute_dft(plan::PlanPtr, X::Ptr{T}, Y::Ptr{T})::Cvoid
    end
end

function unsafe_execute!(p::cWalkPlan{T}, x::StridedArray{T},
                         y::StridedArray{T}) where {T <: fftwComplex}
    walkistride = map(i -> p.istride[i], p.walkdims)
    walkostride = map(i -> p.ostride[i], p.walkdims)
    ax = map(i -> axes(x,i) .- first(axes(x,i)), p.walkdims)
    Elsz, pointerˣ, pointerʸ = sizeof(T), pointer(x), pointer(y)
    for ind in Iterators.product(ax...)
        pointerˣ′ = pointerˣ + Elsz * sum(walkistride .* ind)
        pointerʸ′ = pointerʸ + Elsz * sum(walkostride .* ind)
        unsafe_walk_execute!(p, pointerˣ′, pointerʸ′)
    end
end

function mul!(y::StridedArray{T,N}, p::cWalkPlan{T}, x::StridedArray{T,N}) where {T,N}
    assert_applicable(p, x, y)
    unsafe_execute!(p, x, y)
    y
end

function *(p::cWalkPlan{T,K,inplace}, x::StridedArray{T,N}) where {T,K,N,inplace}
    assert_applicable(p, x)
    y = inplace ? x : Array{T,N}(undef, p.osz)
    unsafe_execute!(p, x, y)
    y
end

const MaybeHighDimension{T} = Union{StridedArray{T,3}, StridedArray{T,4}, StridedArray{T,5}, StridedArray{T,6}, StridedArray{T,7}}
for inplace in (false,true)
    for (f,direction) in ((:fft,FORWARD), (:bfft,BACKWARD))
        plan_f = inplace ? Symbol("plan_",f,"!") : Symbol("plan_",f)
        @eval $plan_f(X::MaybeHighDimension{<:fftwComplex}, region; kws...) = $plan_f(X, Tuple(region); kws...)
        @eval $plan_f(X::MaybeHighDimension{<:fftwComplex}; kws...) = $plan_f(X, ntuple(identity, ndims(X)); kws...)
        @eval $plan_f(X::MaybeHighDimension{<:fftwComplex}, region::Union{Integer,Tuple};
                    flags::Integer=ESTIMATE,
                    timelimit::Real=NO_TIMELIMIT) = begin
            T, N = eltype(X), ndims(X)
            X′ = $inplace ? X : FakeArray{T}(size(X))
            N <= length(region) + 1 && return cFFTWPlan{T,$direction,$inplace,N}(X, X′, region, flags, timelimit)
            cWalkPlan{T,$direction,$inplace}(X, X′, region)
        end
    end
end
