import Base: IndexStyle, getindex, setindex!, eltype, \, similar, copy, real, read!

export PaddedRFFTArray, plan_rfft!, rfft!, plan_irfft!, plan_brfft!, brfft!, irfft!


# This struct reinterprets the `data` array to a Complex or Float array, depending on `eltype(data)`
# It is used internally with the PaddedRFFTArray in place of `Base.ReinterpretArray` 
# ReinterpretArray has some performance issues when reinterprreting a Complex array to Real
struct ComplexOrRealReinterpretArray{T<:fftwNumber,N,A<:DenseArray{<:fftwNumber,N},B<:Ptr} <: DenseArray{T,N}
    data::A # Either a real or complex array
    _unsafe_pointer::B # Pointer to the `data` array, but converted to a different type representation.

    function ComplexOrRealReinterpretArray(rarray::DenseArray{T,N}) where {T<:fftwReal,N}
        ptr = unsafe_convert(Ptr{Complex{T}}, pointer(rarray))
        return new{Complex{T},N,typeof(rarray),typeof(ptr)}(rarray,ptr)
    end

    function ComplexOrRealReinterpretArray(carray::DenseArray{T,N}) where {T<:fftwComplex,N}
        FT = T === ComplexF64 ? Float64 : Float32
        ptr = unsafe_convert(Ptr{FT}, pointer(carray))
        return new{FT,N,typeof(carray),typeof(ptr)}(carray,ptr)
    end
end

const RealReinterpretArray{N} = ComplexOrRealReinterpretArray{<:fftwReal,N,<:DenseArray{<:fftwComplex,N}}
const ComplexReinterpretArray{N} = ComplexOrRealReinterpretArray{<:fftwComplex,N,<:DenseArray{<:fftwReal,N}}

@inline size_convertion(::RealReinterpretArray,i::Integer) = 2i
@inline size_convertion(::ComplexReinterpretArray,i::Integer) = i÷2

IndexStyle(::Type{T}) where {T<:ComplexOrRealReinterpretArray} = IndexLinear()

Base.size(a::ComplexOrRealReinterpretArray) = 
    ntuple(i->(i == 1 ? size_convertion(a,size(a.data)[i]) : size(a.data)[i]),Val(ndims(a.data)))

@inline function getindex(a::ComplexOrRealReinterpretArray,i::Integer)
    data = a.data
    @boundscheck checkbounds(a,i)
    GC.@preserve data r = unsafe_load(a._unsafe_pointer, i)
    return r
end

@inline function setindex!(a::ComplexOrRealReinterpretArray,v,i::Integer)
    data = a.data
    @boundscheck checkbounds(a,i)
    GC.@preserve data unsafe_store!(a._unsafe_pointer,v, i)
    return a
end

Base.unsafe_convert(p::Type{Ptr{T}}, a::ComplexOrRealReinterpretArray{T,N}) where {T,N} = Base.unsafe_convert(p,a.data)

Base.elsize(::Type{<:ComplexOrRealReinterpretArray{T,N}}) where {T,N} = sizeof(T)

complex_or_real_reinterpret(a::AbstractArray) = ComplexOrRealReinterpretArray(a)
complex_or_real_reinterpret(a::ComplexOrRealReinterpretArray) = a.data

# At the time this code was written the new `ReinterpretArray` in Base had some performace issues.
# Those issues were bypassed with the usage of our simplified version of ReinterpretArray above. 
# Hopefully, once the performance issues with ReinterpretArray
# are solved we can just use Base.ReinterpretArray directly.

struct PaddedRFFTArray{T<:fftwReal,N,R,C,L,Nm1} <: DenseArray{Complex{T},N}
    data::R
    r::SubArray{T,N,R,Tuple{Base.OneTo{Int},Vararg{Base.Slice{Base.OneTo{Int}},Nm1}},L} # Real view skipping padding
    c::C

    function PaddedRFFTArray{T,N}(rr::DenseArray{T,N},nx::Int) where {T<:fftwReal,N}
        fsize = size(rr)[1]
        iseven(fsize) || throw(
            ArgumentError("First dimension of allocated array must have even number of elements"))
        (nx == fsize-2 || nx == fsize-1) || throw(
            ArgumentError("Number of elements on the first dimension of array must be either 1 or 2 less than the number of elements on the first dimension of the allocated array"))
        c = complex_or_real_reinterpret(rr)
        r = view(rr, Base.OneTo(nx), ntuple(i->Colon(),Val(N-1))...)
        return  new{T, N, typeof(rr), typeof(c), N===1, N-1}(rr,r,c)
    end # function

    function PaddedRFFTArray{T,N}(c::DenseArray{Complex{T},N},nx::Int) where {T<:fftwReal,N}
        rr = complex_or_real_reinterpret(c)
        fsize = size(rr)[1]
        (nx == fsize-2 || nx == fsize-1) || throw(
            ArgumentError("Number of elements on the first dimension of array must be either 1 or 2 less than the number of elements on the first dimension of the allocated array"))
        r = view(rr, Base.OneTo(nx), ntuple(i->Colon(),Val(N-1))...)
        return  new{T, N, typeof(rr), typeof(c), N===1, N-1}(rr,r,c)
    end # function

end # struct

PaddedRFFTArray(a::DenseArray{<:Union{T,Complex{T}},N},nx::Int) where {T<:fftwReal,N} =
    PaddedRFFTArray{T,N}(a,nx)

function PaddedRFFTArray{T}(ndims::Vararg{Integer,N}) where {T,N}
    fsize = (ndims[1]÷2 + 1)*2
    a = zeros(T,(fsize, ndims[2:end]...))
    PaddedRFFTArray{T,N}(a, ndims[1])
end

PaddedRFFTArray{T}(ndims::NTuple{N,Integer}) where {T,N} =
    PaddedRFFTArray{T}(ndims...)
    
PaddedRFFTArray(ndims::Vararg{Integer,N}) where N = 
    PaddedRFFTArray{Float64}(ndims...)
    
PaddedRFFTArray(ndims::NTuple{N,Integer}) where N = 
    PaddedRFFTArray{Float64}(ndims...)

function PaddedRFFTArray(a::AbstractArray{T,N}) where {T<:fftwReal,N}
    t = PaddedRFFTArray{T}(size(a))
    @inbounds copyto!(t.r, a) 
    return t
end

copy(S::PaddedRFFTArray) = PaddedRFFTArray(copy(S.data),size(S.r,1))

similar(f::PaddedRFFTArray,::Type{T},dims::Tuple{Vararg{Int,N}}) where {T, N} =
    PaddedRFFTArray{T}(dims) 

similar(f::PaddedRFFTArray{T,N,L},dims::NTuple{N2,Int}) where {T,N,L,N2} =
    PaddedRFFTArray{T}(dims) 

similar(f::PaddedRFFTArray,::Type{T}) where {T} =
    PaddedRFFTArray{T}(size(f.r)) 

similar(f::PaddedRFFTArray{T,N}) where {T,N} = 
    PaddedRFFTArray{T,N}(similar(f.data), size(f.r,1)) 

size(S::PaddedRFFTArray) =
    size(S.c)

IndexStyle(::Type{T}) where {T<:PaddedRFFTArray} = 
    IndexLinear()

Base.@propagate_inbounds getindex(A::PaddedRFFTArray,i::Integer) = 
    getindex(A.c,i)

Base.@propagate_inbounds setindex!(A::PaddedRFFTArray,x, i::Integer) = 
    setindex!(A.c,x,i)

Base.unsafe_convert(p::Type{Ptr{Complex{T}}}, a::PaddedRFFTArray{T,N}) where {T,N} = Base.unsafe_convert(p,a.c)

Base.elsize(::Type{<:PaddedRFFTArray{T,N}}) where {T,N} = sizeof(Complex{T})


function PaddedRFFTArray(stream::IO, dims)
    field = PaddedRFFTArray(dims)
    return read!(stream,field)
end

function PaddedRFFTArray{T}(stream::IO, dims) where T
    field = PaddedRFFTArray{T}(dims)
    return read!(stream,field)
end

function read!(file::AbstractString, field::PaddedRFFTArray)
    open(file) do io 
       return read!(io,field) 
    end
end

# Read a binary file of an unpaded array directly to a PaddedRFFT array, without the need
# of the creation of a intermediary Array. If the data is already padded then the user
# should just use PaddedRFFTArray{T}(read("file",unpaddeddim),nx)
function read!(stream::IO, field::PaddedRFFTArray{T,N,L}) where {T,N,L}
    rr = field.data
    dims = size(field.r)
    nx = dims[1]
    nb = sizeof(T)*nx
    npencils = prod(dims)÷nx
    npad = iseven(nx) ? 2 : 1
    for i=0:(npencils-1)
        unsafe_read(stream,Ref(rr,Int((nx+npad)*i+1)),nb)
    end
    return field
end


###########################################################################################
# Foward plans

function plan_rfft!(X::PaddedRFFTArray{T,N}, region;
                   flags::Integer=ESTIMATE,
                   timelimit::Real=NO_TIMELIMIT) where {T<:fftwReal,N}

    (1 in region) || throw(ArgumentError("The first dimension must always be transformed"))
    return rFFTWPlan{T,FORWARD,true,N}(X.r, X.c, region, flags, timelimit)
end

plan_rfft!(f::PaddedRFFTArray;kws...) = plan_rfft!(f, 1:ndims(f); kws...)

*(p::rFFTWPlan{T,FORWARD,true,N},f::PaddedRFFTArray{T,N}) where {T<:fftwReal,N} = 
    (mul!(f.c, p, f.r); f)

rfft!(f::PaddedRFFTArray, region=1:ndims(f)) = plan_rfft!(f, region) * f

function rfft!(r::SubArray{<:fftwReal}, region=1:ndims(r)) 
    f = PaddedRFFTArray(parent(r),size(r,1))
    plan_rfft!(f, region) * f
end

function rfft!(r::DenseArray{<:fftwReal}, region=1:ndims(r)) 
    f = PaddedRFFTArray(r)
    plan_rfft!(f, region) * f
end

function \(p::rFFTWPlan{T,FORWARD,true,N},f::PaddedRFFTArray{T,N}) where {T<:fftwReal,N}
    isdefined(p,:pinv) || (p.pinv = plan_irfft!(f,p.region))
    return p.pinv * f
end


##########################################################################################
# Inverse plans

function plan_brfft!(X::PaddedRFFTArray{T,N}, region;
                    flags::Integer=ESTIMATE,
                    timelimit::Real=NO_TIMELIMIT) where {T<:fftwReal,N}
    (1 in region) || throw(ArgumentError("The first dimension must always be transformed"))
    return rFFTWPlan{Complex{T},BACKWARD,true,N}(X.c, X.r, region, flags,timelimit)
end

plan_brfft!(f::PaddedRFFTArray;kws...) = plan_brfft!(f,1:ndims(f);kws...)

*(p::rFFTWPlan{Complex{T},BACKWARD,true,N},f::PaddedRFFTArray{T,N}) where {T<:fftwReal,N} = 
    (mul!(f.r, p, f.c); f.r)

brfft!(f::PaddedRFFTArray, region=1:ndims(f)) = plan_brfft!(f, region) * f

function brfft!(f::PaddedRFFTArray, i::Integer) 
    if i == size(f.r,1) # Assume `i` is the same as `d` in the brfft!(c::DenseArray{<:fftComplex}, d::Integer, region) defined below
        return brfft!(f,1:ndims(f))
    else # Assume `i` is specifying the region. `plan_brfft!` will throw an error if i != 1
        return brfft!(f,(i,))
    end
end

function brfft!(c::DenseArray{<:fftwComplex}, d::Integer, region=1:ndims(c))
    f = PaddedRFFTArray(c,d)
    plan_brfft!(f, region) * f
end

function plan_irfft!(x::PaddedRFFTArray{T,N}, region; kws...) where {T,N}
    ScaledPlan(plan_brfft!(x, region; kws...),normalization(T, size(x.r), region))
end

plan_irfft!(f::PaddedRFFTArray;kws...) = plan_irfft!(f,1:ndims(f);kws...)

*(p::ScaledPlan{Complex{T},rFFTWPlan{Complex{T},BACKWARD,true,N}},f::PaddedRFFTArray{T,N}) where {T,N} = begin
    p.p * f
    rmul!(f.data, p.scale)
    f.r
end

irfft!(f::PaddedRFFTArray, region=1:ndims(f)) = plan_irfft!(f,region) * f

function irfft!(f::PaddedRFFTArray, i::Integer) 
    if i == size(f.r,1) # Assume `i` is the same as `d` in the irfft!(c::DenseArray{<:fftComplex}, d::Integer, region) defined below
        return irfft!(f,1:ndims(f))
    else # Assume `i` is specifying the region. `plan_irfft!` will throw an error if i != 1
        return irfft!(f,(i,))
    end
end

function irfft!(c::DenseArray{<:fftwComplex}, d::Integer, region=1:ndims(c))
    f = PaddedRFFTArray(c,d)
    plan_irfft!(f, region) * f
end
