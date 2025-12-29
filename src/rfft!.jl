import Base: IndexStyle, getindex, setindex!, eltype, \, similar, copy, real, read!

export PaddedRFFTArray, plan_rfft!, rfft!, plan_irfft!, plan_brfft!, brfft!, irfft!



# As the time this code was written the new `ReinterpretArray` introduced in
# Julia v0.7 had major performace issues. Those issues were bypassed with the usage of the
# custom getindex and setindex! below. Hopefully, once the performance issues with ReinterpretArray
# are solved we can just index the reinterpret array directly.

struct PaddedRFFTArray{T<:fftwReal,N,Nm1,L} <: DenseArray{Complex{T},N}
    data::Array{T,N}
    r::SubArray{T,N,Array{T,N},Tuple{Base.OneTo{Int},Vararg{Base.Slice{Base.OneTo{Int}},Nm1}},L} # Real view skipping padding
    c::Base.ReinterpretArray{Complex{T},N,T,Array{T,N}}

    function PaddedRFFTArray{T,N,Nm1,L}(rr::Array{T,N},nx::Int) where {T<:fftwReal,N,Nm1,L}
        fsize = size(rr)[1]
        iseven(fsize) || throw(
            ArgumentError("First dimension of allocated array must have even number of elements"))
        (nx == fsize-2 || nx == fsize-1) || throw(
            ArgumentError("Number of elements on the first dimension of array must be either 1 or 2 less than the number of elements on the first dimension of the allocated array"))
        c = reinterpret(Complex{T}, rr)
        r = view(rr, Base.OneTo(nx), ntuple(i->Colon(),Val(Nm1))...)
        return  new{T, N, Nm1, L}(rr,r,c)
    end # function
end # struct

@generated function PaddedRFFTArray{T,N}(rr::Array{T,N},nx::Int) where {T<:fftwReal,N}
    :(PaddedRFFTArray{T,N,$(N-1),$(N === 1 ? true : false)}(rr,nx))
end

@inline real(S::PaddedRFFTArray) = S.r

@inline complex_view(S::PaddedRFFTArray) = S.c

@inline data(S::PaddedRFFTArray) = S.data

copy(S::PaddedRFFTArray) = PaddedRFFTArray(copy(data(S)),size(real(S),1))

similar(f::PaddedRFFTArray,::Type{T},dims::Tuple{Vararg{Int,N}}) where {T, N} =
    PaddedRFFTArray{T}(dims) 
similar(f::PaddedRFFTArray{T,N,L},dims::NTuple{N2,Int}) where {T,N,L,N2} =
    PaddedRFFTArray{T}(dims) 
similar(f::PaddedRFFTArray,::Type{T}) where {T} =
    PaddedRFFTArray{T}(size(real(f))) 
similar(f::PaddedRFFTArray{T,N}) where {T,N} = 
    PaddedRFFTArray{T,N}(similar(data(f)), size(real(f),1)) 

size(S::PaddedRFFTArray) =
    size(complex_view(S))

IndexStyle(::Type{T}) where {T<:PaddedRFFTArray} = 
    IndexLinear()

@inline function getindex(A::PaddedRFFTArray{T,N}, i2::Integer) where {T,N}
    d = data(A)
    i = 2i2
    @boundscheck checkbounds(d,i)
    @inbounds begin 
        return Complex{T}(d[i-1],d[i])
    end
end    

@inline @generated function getindex(A::PaddedRFFTArray{T,N}, I2::Vararg{Integer,N}) where {T,N}
    ip = :(2*I2[1])
    t = Expr(:tuple)
    for i=2:N
        push!(t.args,:(I2[$i]))
    end
    quote
        d = data(A)
        i = $ip
        I = $t
        @boundscheck checkbounds(d,i,I...)
        @inbounds begin 
            return Complex{T}(d[i-1,I...],d[i,I...])
        end
    end
end

@inline function setindex!(A::PaddedRFFTArray{T,N},x, i2::Integer) where {T,N}
    d = data(A)
    i = 2i2
    @boundscheck checkbounds(d,i)
    @inbounds begin 
        d[i-1] = real(x)
        d[i] = imag(x)
    end
    A
end

@inline @generated function setindex!(A::PaddedRFFTArray{T,N}, x, I2::Vararg{Integer,N}) where {T,N}
    ip = :(2*I2[1])
    t = Expr(:tuple)
    for i=2:N
        push!(t.args,:(I2[$i]))
    end
    quote
        d = data(A)
        i = $ip
        I = $t
        @boundscheck checkbounds(d,i,I...)
        @inbounds begin 
            d[i-1,I...] = real(x)
            d[i,I...] = imag(x)
        end
        A
    end
end

PaddedRFFTArray(rr::Array{T,N},nx::Int) where {T<:fftwReal,N} = PaddedRFFTArray{T,N}(rr,nx)

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

function PaddedRFFTArray{T}(a::AbstractArray{<:Real,N}) where {T<:fftwReal,N}
    t = PaddedRFFTArray{T}(size(a))
    @inbounds copyto!(t.r, a) 
    return t
end

PaddedRFFTArray(a::AbstractArray{<:Real}) = PaddedRFFTArray{Float64}(a)

function PaddedRFFTArray(stream, dims)
    field = PaddedRFFTArray(dims)
    return read!(stream,field)
end

function PaddedRFFTArray{T}(stream, dims) where T
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
    rr = data(field)
    dims = size(real(field))
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
    return rFFTWPlan{T,FORWARD,true,N}(real(X), complex_view(X), region, flags, timelimit)
end

plan_rfft!(f::PaddedRFFTArray;kws...) = plan_rfft!(f, 1:ndims(f); kws...)

*(p::rFFTWPlan{T,FORWARD,true,N},f::PaddedRFFTArray{T,N}) where {T<:fftwReal,N} = 
    (mul!(complex_view(f), p, real(f)); f)

rfft!(f::PaddedRFFTArray, region=1:ndims(f)) = plan_rfft!(f, region) * f

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
    return rFFTWPlan{Complex{T},BACKWARD,true,N}(complex_view(X), real(X), region, flags,timelimit)
end

plan_brfft!(f::PaddedRFFTArray;kws...) = plan_brfft!(f,1:ndims(f);kws...)

*(p::rFFTWPlan{Complex{T},BACKWARD,true,N},f::PaddedRFFTArray{T,N}) where {T<:fftwReal,N} = 
    (mul!(real(f), p, complex_view(f)); real(f))

brfft!(f::PaddedRFFTArray, region=1:ndims(f)) = plan_brfft!(f, region) * f

function plan_irfft!(x::PaddedRFFTArray{T,N}, region; kws...) where {T,N}
    ScaledPlan(plan_brfft!(x, region; kws...),normalization(T, size(real(x)), region))
end

plan_irfft!(f::PaddedRFFTArray;kws...) = plan_irfft!(f,1:ndims(f);kws...)

*(p::ScaledPlan,f::PaddedRFFTArray) = begin
    p.p * f
    rmul!(data(f), p.scale)
    real(f)
end

irfft!(f::PaddedRFFTArray, region=1:ndims(f)) = plan_irfft!(f,region) * f
