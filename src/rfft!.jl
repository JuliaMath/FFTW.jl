import Base: IndexStyle, getindex, setindex!, eltype, \, similar, copy, real, read!

#For compatibility between Julia v0.6 and v0.7 - begin
if VERSION >= v"0.7-"
    using Base.@gc_preserve
else
    macro gc_preserve(s::Symbol,ex::Expr)
        return esc(ex)
    end
end
#For compatibility between Julia v0.6 and v0.7 - end

export PaddedRFFTArray, plan_rfft!, rfft!, plan_irfft!, plan_brfft!, brfft!, irfft!



# As the time this code was written the new `ReinterpretArray` introduced in
# Julia v0.7 had major performace issues. Those issues were bypassed with the usage of the
# unsafe_wrap for the complex view of the data. As the name sugest, this is
# unsafe: whenever this complex view is called one must be careful to gc_preserve
# the "parent" PaddedRFFTArray. Therefore this view should not be "exported" to the
# user and the PaddedRFFTArray itself behaves like the complex view.
# Since it is not possible to not export a particular field of a exported type,
# a "hack" was used to name the unsafe field "#c", a fieldname that a non-advanced user
# will likely not be able to call.

const c = Symbol("#c") 
@eval struct PaddedRFFTArray{T<:fftwReal,N,L} <: DenseArray{Complex{T},N}
    r::SubArray{T,N,Array{T,N},NTuple{N,UnitRange{Int}},L} # Real view skipping padding
    ($c)::Array{Complex{T},N}

    function PaddedRFFTArray{T,N}(rr::Array{T,N},nx::Int) where {T<:fftwReal,N}
        rrsize = size(rr)
        fsize = rrsize[1]
        iseven(fsize) || throw(
            ArgumentError("First dimension of allocated array must have even number of elements"))
        (nx == fsize-2 || nx == fsize-1) || throw(
            ArgumentError("Number of elements on the first dimension of array must be either 1 or 2 less than the number of elements on the first dimension of the allocated array"))
        fsize = fsize÷2
        csize = (fsize, rrsize[2:end]...)
        if VERSION >= v"0.7-" 
            @gc_preserve rr c = unsafe_wrap(Array{Complex{T},N}, 
                                           reinterpret(Ptr{Complex{T}},pointer(rr)), 
                                           csize)
        else 
            c = reinterpret(Complex{T}, rr, csize)
        end
        rsize = (nx,rrsize[2:end]...)
        r = view(rr,(1:l for l in rsize)...)
        return  @gc_preserve rr new{T, N, N === 1 ? true : false}(r,c)
    end # function
end # struct

@inline real(S::PaddedRFFTArray) = S.r

@inline unsafe_complex_view(S::PaddedRFFTArray) = getfield(S,c)

copy(S::PaddedRFFTArray) = PaddedRFFTArray(copy(parent(real(S))),size(real(S))[1])

similar(f::PaddedRFFTArray,::Type{T},dims::Tuple{Vararg{Int,N}}) where {T, N} =
    PaddedRFFTArray{T}(dims) 
similar(f::PaddedRFFTArray{T,N,L},dims::NTuple{N2,Int}) where {T,N,L,N2} =
    PaddedRFFTArray{T}(dims) 
similar(f::PaddedRFFTArray,::Type{T}) where {T} =
    PaddedRFFTArray{T}(size(real(f))) 
similar(f::PaddedRFFTArray{T,N}) where {T,N} = 
    PaddedRFFTArray{T,N}(similar(parent(real(f))), size(real(f))[1]) 

size(S::PaddedRFFTArray) =
    @gc_preserve S size(unsafe_complex_view(S))

IndexStyle(::Type{T}) where {T<:PaddedRFFTArray} = 
    IndexLinear()

Base.@propagate_inbounds @inline getindex(S::PaddedRFFTArray, i::Int) =
    @gc_preserve S getindex(unsafe_complex_view(S),i)

Base.@propagate_inbounds @inline getindex(S::PaddedRFFTArray{T,N}, I::Vararg{Int, N}) where {T,N} =
    @gc_preserve S getindex(unsafe_complex_view(S),I...)

Base.@propagate_inbounds @inline setindex!(S::PaddedRFFTArray,v,i::Int) = 
    @gc_preserve S setindex!(unsafe_complex_view(S),v,i)

Base.@propagate_inbounds @inline setindex!(S::PaddedRFFTArray{T,N},v,I::Vararg{Int,N}) where {T,N} =
    @gc_preserve S setindex!(unsafe_complex_view(S),v,I...)


PaddedRFFTArray(rr::Array{T,N},nx::Int) where {T<:fftwReal,N} = PaddedRFFTArray{T,N}(rr,nx)

function PaddedRFFTArray{T}(ndims::Vararg{Integer,N}) where {T,N}
    fsize = (ndims[1]÷2 + 1)*2
    a = Array{T,N}((fsize, ndims[2:end]...))
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
    @inbounds copy!(t.r, a) 
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
    rr = parent(field.r)
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
    if flags&ESTIMATE != 0
        @gc_preserve X p = 
            rFFTWPlan{T,FORWARD,true,N}(real(X), unsafe_complex_view(X), region, flags, timelimit)
    else
        x = similar(X)
        @gc_preserve x p =
            rFFTWPlan{T,FORWARD,true,N}(real(x), unsafe_complex_view(x), region, flags, timelimit)
    end
    return p
end

plan_rfft!(f::PaddedRFFTArray;kws...) = plan_rfft!(f, 1:ndims(f); kws...)

*(p::rFFTWPlan{T,FORWARD,true,N},f::PaddedRFFTArray{T,N}) where {T<:fftwReal,N} = 
    (@gc_preserve f A_mul_B!(unsafe_complex_view(f), p, real(f)); f)

rfft!(f::PaddedRFFTArray, region=1:ndims(f)) = plan_rfft!(f, region) * f

function \(p::rFFTWPlan{T,FORWARD,true,N},f::PaddedRFFTArray{T,N}) where {T<:fftwReal,N}
    isdefined(p,:pinv) || (p.pinv = plan_irfft!(f,p.region))
    return p.pinv * f
end


##########################################################################################
# Inverse plans

function plan_brfft!(X::PaddedRFFTArray{T,N}, region;
                    flags::Integer=PRESERVE_INPUT,
                    timelimit::Real=NO_TIMELIMIT) where {T<:fftwReal,N}
    (1 in region) || throw(ArgumentError("The first dimension must always be transformed"))
    if flags&PRESERVE_INPUT != 0
        a = similar(X)
        return @gc_preserve a rFFTWPlan{Complex{T},BACKWARD,true,N}(unsafe_complex_view(a), real(a), region, flags,timelimit)
    else
        return @gc_preserve X rFFTWPlan{Complex{T},BACKWARD,true,N}(unsafe_complex_view(X), real(X), region, flags,timelimit)
    end
end

plan_brfft!(f::PaddedRFFTArray;kws...) = plan_brfft!(f,1:ndims(f);kws...)

*(p::rFFTWPlan{Complex{T},BACKWARD,true,N},f::PaddedRFFTArray{T,N}) where {T<:fftwReal,N} = 
    (@gc_preserve f A_mul_B!(real(f), p, unsafe_complex_view(f)); real(f))

brfft!(f::PaddedRFFTArray, region=1:ndims(f)) = plan_brfft!(f, region) * f

function plan_irfft!(x::PaddedRFFTArray{T,N}, region; kws...) where {T,N}
    ScaledPlan(plan_brfft!(x, region; kws...),normalization(T, size(real(x)), region))
end

plan_irfft!(f::PaddedRFFTArray;kws...) = plan_irfft!(f,1:ndims(f);kws...)

*(p::ScaledPlan,f::PaddedRFFTArray) = begin
    p.p * f
    scale!(parent(real(f)), p.scale)
    real(f)
end

irfft!(f::PaddedRFFTArray, region=1:ndims(f)) = plan_irfft!(f,region) * f




