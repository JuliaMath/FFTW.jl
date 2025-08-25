# This file was formerly a part of Julia. License is MIT: https://julialang.org/license

import Base: show, *, convert, unsafe_convert, size, strides, ndims, pointer
import LinearAlgebra: mul!

"""
    r2r(A, kind [, dims])

Performs a multidimensional real-input/real-output (r2r) transform
of type `kind` of the array `A`, as defined in the FFTW manual.
`kind` specifies either a discrete cosine transform of various types
(`FFTW.REDFT00`, `FFTW.REDFT01`, `FFTW.REDFT10`, or
`FFTW.REDFT11`), a discrete sine transform of various types
(`FFTW.RODFT00`, `FFTW.RODFT01`, `FFTW.RODFT10`, or
`FFTW.RODFT11`), a real-input DFT with halfcomplex-format output
(`FFTW.R2HC` and its inverse `FFTW.HC2R`), or a discrete
Hartley transform (`FFTW.DHT`).  The `kind` argument may be
an array or tuple in order to specify different transform types
along the different dimensions of `A`; `kind[end]` is used
for any unspecified dimensions.  See the FFTW manual for precise
definitions of these transform types, at http://www.fftw.org/doc.

The optional `dims` argument specifies an iterable subset of
dimensions (e.g. an integer, range, tuple, or array) to transform
along. `kind[i]` is then the transform type for `dims[i]`,
with `kind[end]` being used for `i > length(kind)`.

See also [`plan_r2r`](@ref) to pre-plan optimized r2r transforms.
"""
function r2r end

"""
    r2r!(A, kind [, dims])

Same as [`r2r`](@ref), but operates in-place on `A`, which must be
an array of real or complex floating-point numbers.
"""
function r2r! end

"""
    plan_r2r!(A, kind [, dims [, flags [, timelimit [, num_threads]]]])

Similar to [`plan_fft`](@extref `AbstractFFTs.plan_fft`), but corresponds to [`r2r!`](@ref).
"""
function plan_r2r! end

"""
    plan_r2r(A, kind [, dims [, flags [, timelimit [, num_threads]]]])

Pre-plan an optimized r2r transform, similar to [`plan_fft`](@extref `AbstractFFTs.plan_fft`) except that the
transforms (and the first three arguments) correspond to [`r2r`](@ref) and [`r2r!`](@ref), respectively.
"""
function plan_r2r end

## FFT: Implement fft by calling fftw.

# TODO: this is dangerous since it captures runtime data from the compile machine
# n.b. write this as a function to avoid Julia bugs with the runtime cglobal implementation
get_version() = VersionNumber(split(unsafe_string(cglobal(
    (:fftw_version,libfftw3_no_init), UInt8)), ['-', ' '])[2])
const version = get_version()

## Direction of FFT

const FORWARD = -1
const BACKWARD = 1

## FFTW Flags from fftw3.h

const MEASURE         = UInt32(0)
const DESTROY_INPUT   = UInt32(1 << 0)
const UNALIGNED       = UInt32(1 << 1)
const CONSERVE_MEMORY = UInt32(1 << 2)
const EXHAUSTIVE      = UInt32(1 << 3)   # NO_EXHAUSTIVE is default
const PRESERVE_INPUT  = UInt32(1 << 4)   # cancels DESTROY_INPUT
const PATIENT         = UInt32(1 << 5)   # IMPATIENT is default
const ESTIMATE        = UInt32(1 << 6)
const WISDOM_ONLY     = UInt32(1 << 21)
const NO_SIMD = UInt32(1 << 17) # disable SIMD, useful for benchmarking

## R2R transform kinds

const R2HC    = 0
const HC2R    = 1
const DHT     = 2
const REDFT00 = 3
const REDFT01 = 4
const REDFT10 = 5
const REDFT11 = 6
const RODFT00 = 7
const RODFT01 = 8
const RODFT10 = 9
const RODFT11 = 10

let k2s = Dict(R2HC => "R2HC", HC2R => "HC2R", DHT => "DHT", REDFT00 => "REDFT00",
               REDFT01 => "REDFT01", REDFT10 => "REDFT10", REDFT11 => "REDFT11",
               RODFT00 => "RODFT00", RODFT01 => "RODFT01", RODFT10 => "RODFT10",
               RODFT11 => "RODFT11")
    global kind2string
    kind2string(k::Integer) = k2s[Int(k)]
end

# FFTW floating-point types:

const fftwNumber = Union{Float64,Float32,Complex{Float64},Complex{Float32}}
const fftwReal = Union{Float64,Float32}
const fftwComplex = Union{Complex{Float64},Complex{Float32}}
const fftwDouble = Union{Float64,Complex{Float64}}
const fftwSingle = Union{Float32,Complex{Float32}}
const fftwTypeDouble = Union{Type{Float64},Type{Complex{Float64}}}
const fftwTypeSingle = Union{Type{Float32},Type{Complex{Float32}}}

# For ESTIMATE plans, FFTW allows one to pass NULL for the array pointer,
# since it is not written to.  Hence, it is convenient to create an
# array-like type that carries a size and a stride like a "real" array
# but which is converted to C_NULL as a pointer.
struct FakeArray{T,N} <: DenseArray{T,N}
    sz::NTuple{N,Int}
    st::NTuple{N,Int}
end
size(a::FakeArray) = a.sz
strides(a::FakeArray) = a.st
unsafe_convert(::Type{Ptr{T}}, a::FakeArray{T}) where {T} = convert(Ptr{T}, C_NULL)
pointer(a::FakeArray{T}) where {T} = convert(Ptr{T}, C_NULL)
FakeArray{T}(sz::NTuple{N,Int}) where {T,N} = FakeArray{T,N}(sz, colmajorstrides(sz))
FakeArray{T}(sz::Int...) where {T} = FakeArray{T}(sz)
fakesimilar(flags, X, T) = flags & ESTIMATE != 0 ? FakeArray{T}(size(X)) : Array{T}(undef, size(X))
alignment_of(A::FakeArray) = Int32(0)

## Julia wrappers around FFTW functions

# Wisdom

# Import and export wisdom to/from a single file for all precisions,
# which is more user-friendly than requiring the user to call a
# separate routine depending on the fp precision of the plans.  This
# requires a bit of trickiness since we have to (a) use the libc file
# I/O routines with fftw_export_wisdom_to_file/import_wisdom_from_file
# (b) we need 256 bytes of space padding between the wisdoms to work
# around FFTW's internal file i/o buffering [see the BUFSZ constant in
# FFTW's api/import-wisdom-from-file.c file].

@exclusive function export_wisdom(fname::AbstractString)
    f = ccall(:fopen, Ptr{Cvoid}, (Cstring,Cstring), fname, :w)
    systemerror("could not open wisdom file $fname for writing", f == C_NULL)
    ccall((:fftw_export_wisdom_to_file,libfftw3), Cvoid, (Ptr{Cvoid},), f)
    ccall(:fputs, Int32, (Ptr{UInt8},Ptr{Cvoid}), " "^256, f) # no NUL, hence no Cstring
    ccall((:fftwf_export_wisdom_to_file,libfftw3f), Cvoid, (Ptr{Cvoid},), f)
    ccall(:fclose, Cvoid, (Ptr{Cvoid},), f)
end

@exclusive function import_wisdom(fname::AbstractString)
    f = ccall(:fopen, Ptr{Cvoid}, (Cstring,Cstring), fname, :r)
    systemerror("could not open wisdom file $fname for reading", f == C_NULL)
    if ccall((:fftw_import_wisdom_from_file,libfftw3),Int32,(Ptr{Cvoid},),f)==0||
       ccall((:fftwf_import_wisdom_from_file,libfftw3f),Int32,(Ptr{Cvoid},),f)==0
        error("failed to import wisdom from $fname")
    end
    ccall(:fclose, Cvoid, (Ptr{Cvoid},), f)
end

@exclusive function import_system_wisdom()
    if ccall((:fftw_import_system_wisdom,libfftw3), Int32, ()) == 0 ||
       ccall((:fftwf_import_system_wisdom,libfftw3f), Int32, ()) == 0
        error("failed to import system wisdom")
    end
end

@exclusive function forget_wisdom()
    ccall((:fftw_forget_wisdom,libfftw3), Cvoid, ())
    ccall((:fftwf_forget_wisdom,libfftw3f), Cvoid, ())
end

# Threads

# Must only be called after acquiring fftwlock
function _set_num_threads(num_threads::Integer)
    @static if fftw_provider == "mkl"
        _last_num_threads[] = num_threads
    end
    ccall((:fftw_plan_with_nthreads,libfftw3), Cvoid, (Int32,), num_threads)
    ccall((:fftwf_plan_with_nthreads,libfftw3f), Cvoid, (Int32,), num_threads)
end

@exclusive set_num_threads(num_threads::Integer) = _set_num_threads(num_threads)

function get_num_threads()
    @static if fftw_provider == "fftw"
        ccall((:fftw_planner_nthreads,libfftw3), Cint, ())
    else
        _last_num_threads[]
    end
end

@exclusive function set_num_threads(f::Function, num_threads::Integer)
    orig_num_threads = get_num_threads()
    _set_num_threads(num_threads)
    try
        f()
    finally
        _set_num_threads(orig_num_threads)
    end
end

# pointer type for fftw_plan (opaque pointer)

struct fftw_plan_struct end
const PlanPtr = Ptr{fftw_plan_struct}

# Planner timelimits

const NO_TIMELIMIT = -1.0 # from fftw3.h

# only call these when fftwlock is held:
unsafe_set_timelimit(precision::fftwTypeDouble,seconds) =
    ccall((:fftw_set_timelimit,libfftw3), Cvoid, (Float64,), seconds)
unsafe_set_timelimit(precision::fftwTypeSingle,seconds) =
    ccall((:fftwf_set_timelimit,libfftw3f), Cvoid, (Float64,), seconds)
@exclusive set_timelimit(precision, seconds) = unsafe_set_timelimit(precision, seconds)

# Array alignment mod 16:
#   FFTW plans may depend on the alignment of the array mod 16 bytes,
#   i.e. the address mod 16 of the first element of the array, in order
#   to exploit SIMD operations.  Julia arrays are, by default, aligned
#   to 16-byte boundaries (address mod 16 == 0), but this may not be
#   true for data imported from external C code, or for SubArrays.
#   Use the undocumented routine fftw_alignment_of to determine the
#   alignment of a given pointer modulo whatever FFTW needs; this
#   function will be documented in FFTW 3.3.4.


@static if fftw_provider == "mkl"
    alignment_of(A::StridedArray{<:fftwDouble}) =
        convert(Int32, convert(Int64, pointer(A)) % 16)
    alignment_of(A::StridedArray{<:fftwSingle}) =
        convert(Int32, convert(Int64, pointer(A)) % 16)
else
    alignment_of(A::StridedArray{T}) where {T<:fftwDouble} =
        ccall((:fftw_alignment_of, libfftw3), Int32, (Ptr{T},), A)
    alignment_of(A::StridedArray{T}) where {T<:fftwSingle} =
        ccall((:fftwf_alignment_of, libfftw3f), Int32, (Ptr{T},), A)
end

# FFTWPlan (low-level)

# low-level storage of the FFTW plan, along with the information
# needed to determine whether it is applicable.   We need to put
# this into a type to support a finalizer on the fftw_plan.
# K is FORWARD/BACKWARD for forward/backward or r2c/c2r plans, respectively.
# For r2r plans, the field kinds::K is a tuple/vector of the transform kinds along each dimension.
abstract type FFTWPlan{T<:fftwNumber,K,inplace} <: Plan{T} end
for P in (:cFFTWPlan, :rFFTWPlan) # complex, r2c/c2r
    @eval begin
        mutable struct $P{T<:fftwNumber,K,inplace,N,G} <: FFTWPlan{T,K,inplace}
            plan::PlanPtr
            sz::NTuple{N,Int} # size of array on which plan operates (Int tuple)
            osz::NTuple{N,Int} # size of output array (Int tuple)
            istride::NTuple{N,Int} # strides of input
            ostride::NTuple{N,Int} # strides of output
            ialign::Int32 # alignment mod 16 of input
            oalign::Int32 # alignment mod 16 of input
            flags::UInt32 # planner flags
            region::G # region (iterable) of dims that are transformed
            pinv::ScaledPlan
            function $P{T,K,inplace,N,G}(plan::PlanPtr, flags::Integer, R::G,
                                         X::StridedArray{T,N}, Y::StridedArray) where {T<:fftwNumber,K,inplace,N,G}
                p = new(plan, size(X), size(Y), strides(X), strides(Y),
                        alignment_of(X), alignment_of(Y), flags, R)
                finalizer(maybe_destroy_plan, p)
                p
            end
        end

        function $P{T,K,inplace,N}(plan::PlanPtr, flags::Integer, R::G,
                                   X::StridedArray{T,N},
                                   Y::StridedArray) where {T<:fftwNumber,K,inplace,N,G}
            $P{T,K,inplace,N,G}(plan, flags, R, X, Y)
        end
    end
end

mutable struct r2rFFTWPlan{T<:fftwNumber,K,inplace,N,G} <: FFTWPlan{T,K,inplace}
    plan::PlanPtr
    sz::NTuple{N,Int} # size of array on which plan operates (Int tuple)
    osz::NTuple{N,Int} # size of output array (Int tuple)
    istride::NTuple{N,Int} # strides of input
    ostride::NTuple{N,Int} # strides of output
    ialign::Int32 # alignment mod 16 of input
    oalign::Int32 # alignment mod 16 of input
    flags::UInt32 # planner flags
    region::G # region (iterable) of dims that are transformed
    kinds::K
    pinv::ScaledPlan
    function r2rFFTWPlan{T,K,inplace,N,G}(plan::PlanPtr, flags::Integer, R::G,
                                 X::StridedArray{T,N},
                                 Y::StridedArray, kinds::K) where {T<:fftwNumber,K,inplace,N,G}
        p = new(plan, size(X), size(Y), strides(X), strides(Y),
                alignment_of(X), alignment_of(Y), flags, R, kinds)
        finalizer(maybe_destroy_plan, p)
        p
    end
end

function r2rFFTWPlan{T,K,inplace,N}(plan::PlanPtr, flags::Integer, R::G,
                           X::StridedArray{T,N},
                           Y::StridedArray, kinds::K) where {T<:fftwNumber,K,inplace,N,G}
    r2rFFTWPlan{T,K,inplace,N,G}(plan, flags, R, X, Y, kinds)
end

size(p::FFTWPlan) = p.sz

unsafe_convert(::Type{PlanPtr}, p::FFTWPlan) = p.plan

#################################################################################################
# We have to be careful about destroying plans in an FFTWPlan finalizer, because the
# garbage-collector may run at unexpected times.  For example, it may run during the spawnloop
# callback while the FFTW planner is executing.  Under these circumstances, the garbage collector
# may deadlock if it needs to acquire the fftwlock (already held by the planning thread via @exclusive).
# Instead, if islocked(fftwlock), we defer plan destruction until the fftwlock is released, by
# pushing the plan to be destroyed to the deferred_destroy_plans (which itself is protected by a lock).
# This is accomplished by the maybe_destroy_plan function, which is used as the plan finalizer.

# these functions should only be called while the fftwlock is held
unsafe_destroy_plan(@nospecialize(plan::FFTWPlan{<:fftwDouble})) =
    ccall((:fftw_destroy_plan,libfftw3), Cvoid, (PlanPtr,), plan)
unsafe_destroy_plan(@nospecialize(plan::FFTWPlan{<:fftwSingle})) =
    ccall((:fftwf_destroy_plan,libfftw3f), Cvoid, (PlanPtr,), plan)

const deferred_destroy_lock = ReentrantLock() # lock protecting the deferred_destroy_plans list
const deferred_destroy_plans = FFTWPlan[]

function destroy_deferred()
    lock(deferred_destroy_lock)
    try
        # need trylock here to avoid potential deadlocks if another
        # @exclusive function has just grabbed the lock; in that case,
        # we'll do nothing (the other function will eventually run destroy_deferred).
        if !isempty(deferred_destroy_plans) && trylock(fftwlock)
            try
                @inline foreach(unsafe_destroy_plan, deferred_destroy_plans)
                empty!(deferred_destroy_plans)
            finally
                unlock(fftwlock)
            end
        end
    finally
        unlock(deferred_destroy_lock)
    end
end

@exclusive destroy_plan(plan::FFTWPlan) = unsafe_destroy_plan(plan)

function maybe_destroy_plan(plan::FFTWPlan)
    # need to acquire deferred_destroy_lock before trylock to avoid a memory leak
    # (race to avoid: trylock == false, other thread unlocks and calls destroy_deferred,
    #                 then we push a deferred plan that may never get freed)
    # Also, we need to use a trylock spinloop rather than lock(deferred_destroy_lock),
    # since task switches aren't permitted in finalizers.   This has suboptimal efficiency,
    # but we shouldn't waste too many cycles since destroying plans is quick and contention
    # should be rare.
    while !trylock(deferred_destroy_lock)
        # Need a safepoint in here because without it, this loop blocks forward progress in the GC and can deadlock
        # the program.
        GC.safepoint()
    end
    try
        # note: fftwlock is re-entrant, so trylock will succeed here if we
        # are in the task that holds the planner lock.  That's okay — 
        # re-entrant calls to destroy_plan in the planner thread
        # should be fine as long as we didn't call make_planner_thread_safe.
        if trylock(fftwlock)
            try
                unsafe_destroy_plan(plan)
            finally
                unlock(fftwlock)
            end
        else
            push!(deferred_destroy_plans, plan)
        end
    finally
        unlock(deferred_destroy_lock)
    end
end

#################################################################################################

cost(plan::FFTWPlan{<:fftwDouble}) =
    ccall((:fftw_cost,libfftw3), Float64, (PlanPtr,), plan)
cost(plan::FFTWPlan{<:fftwSingle}) =
    ccall((:fftwf_cost,libfftw3f), Float64, (PlanPtr,), plan)

@exclusive function arithmetic_ops(plan::FFTWPlan{<:fftwDouble})
    add, mul, fma = Ref(0.0), Ref(0.0), Ref(0.0)
    ccall((:fftw_flops,libfftw3), Cvoid,
          (PlanPtr,Ref{Float64},Ref{Float64},Ref{Float64}), plan, add, mul, fma)
    return (round(Int64, add[]), round(Int64, mul[]), round(Int64, fma[]))
end
@exclusive function arithmetic_ops(plan::FFTWPlan{<:fftwSingle})
    add, mul, fma = Ref(0.0), Ref(0.0), Ref(0.0)
    ccall((:fftwf_flops,libfftw3f), Cvoid,
          (PlanPtr,Ref{Float64},Ref{Float64},Ref{Float64}), plan, add, mul, fma)
    return (round(Int64, add[]), round(Int64, mul[]), round(Int64, fma[]))
end
flops(plan::FFTWPlan) = let ops = arithmetic_ops(plan)
    ops[1] + ops[2] + 2 * ops[3] # add + mul + 2*fma
end

# Pretty-printing plans

function showfftdims(io, sz::Dims, istride::Dims, T)
    if isempty(sz)
        print(io, "0-dimensional")
    elseif length(sz) == 1
        print(io, sz[1], "-element")
    else
        print(io, join(sz, "×"))
    end
    if istride == colmajorstrides(sz)
        print(io, " array of ", T)
    else
        print(io, " $istride-strided array of ", T)
    end
end

# The sprint_plan function was released in FFTW 3.3.4, but MKL versions
# claiming to be FFTW 3.3.4 still don't seem to have this function.
const has_sprint_plan = version >= v"3.3.4" && fftw_provider == "fftw"

@static if has_sprint_plan
    sprint_plan_(plan::FFTWPlan{<:fftwDouble}) =
        ccall((:fftw_sprint_plan,libfftw3), Ptr{UInt8}, (PlanPtr,), plan)
    sprint_plan_(plan::FFTWPlan{<:fftwSingle}) =
        ccall((:fftwf_sprint_plan,libfftw3f), Ptr{UInt8}, (PlanPtr,), plan)
    function sprint_plan(plan::FFTWPlan)
        p = sprint_plan_(plan)
        str = unsafe_string(p)
        Libc.free(p)
        return str
    end
else
    sprint_plan(plan::FFTWPlan) = ""
end

function show(io::IO, p::cFFTWPlan{T,K,inplace}) where {T,K,inplace}
    print(io, inplace ? "FFTW in-place " : "FFTW ",
          K < 0 ? "forward" : "backward", " plan for ")
    showfftdims(io, p.sz, p.istride, T)
    has_sprint_plan && print(io, "\n", sprint_plan(p))
end

function show(io::IO, p::rFFTWPlan{T,K,inplace}) where {T,K,inplace}
    print(io, inplace ? "FFTW in-place " : "FFTW ",
          K < 0 ? "real-to-complex" : "complex-to-real",
          " plan for ")
    showfftdims(io, p.sz, p.istride, T)
    has_sprint_plan && print(io, "\n", sprint_plan(p))
end

function show(io::IO, p::r2rFFTWPlan{T,K,inplace}) where {T,K,inplace}
    kinds = p.kinds
    print(io, inplace ? "FFTW in-place r2r " : "FFTW r2r ")
    if isempty(kinds)
        print(io, "0-dimensional")
    elseif kinds == ntuple(i -> kinds[1], length(kinds))
        print(io, kind2string(kinds[1]))
        if length(kinds) > 1
            print(io, "^", length(kinds))
        end
    else
        print(io, join(map(kind2string, kinds), "×"))
    end
    print(io, " plan for ")
    showfftdims(io, p.sz, p.istride, T)
    has_sprint_plan && print(io, "\n", sprint_plan(p))
end

# Check whether a FFTWPlan is applicable to a given input array, and
# throw an informative error if not:
function assert_applicable(p::FFTWPlan{T}, X::StridedArray{T}) where T
    if size(X) != p.sz
        throw(ArgumentError("FFTW plan applied to wrong-size array"))
    elseif strides(X) != p.istride
        throw(ArgumentError("FFTW plan applied to wrong-strides array"))
    elseif alignment_of(X) != p.ialign && p.flags & UNALIGNED == 0
        throw(ArgumentError("FFTW plan applied to array with wrong memory alignment"))
    end
end

function assert_applicable(p::FFTWPlan{T,K,inplace}, X::StridedArray{T}, Y::StridedArray) where {T,K,inplace}
    assert_applicable(p, X)
    if size(Y) != p.osz
        throw(ArgumentError("FFTW plan applied to wrong-size output"))
    elseif strides(Y) != p.ostride
        throw(ArgumentError("FFTW plan applied to wrong-strides output"))
    elseif alignment_of(Y) != p.oalign && p.flags & UNALIGNED == 0
        throw(ArgumentError("FFTW plan applied to output with wrong memory alignment"))
    elseif inplace != (pointer(X) == pointer(Y))
        throw(ArgumentError(string("FFTW ",
                                   inplace ? "in-place" : "out-of-place",
                                   " plan applied to ",
                                   inplace ? "out-of-place" : "in-place",
                                   " data")))
    end
end

# strides for a column-major (Julia-style) array of size == sz
colmajorstrides(::Tuple{}) = ()
colmajorstrides(sz) = _colmajorstrides(1, sz...)
@inline _colmajorstrides(p, sz1, tail...) = (p, _colmajorstrides(p*sz1, tail...)...)
_colmajorstrides(p) = ()

# Execute

unsafe_execute!(plan::FFTWPlan{<:fftwDouble}) =
    ccall((:fftw_execute,libfftw3), Cvoid, (PlanPtr,), plan)

unsafe_execute!(plan::FFTWPlan{<:fftwSingle}) =
    ccall((:fftwf_execute,libfftw3f), Cvoid, (PlanPtr,), plan)

unsafe_execute!(plan::cFFTWPlan{T},
                X::StridedArray{T}, Y::StridedArray{T}) where {T<:fftwDouble} =
    ccall((:fftw_execute_dft,libfftw3), Cvoid,
          (PlanPtr,Ptr{T},Ptr{T}), plan, X, Y)

unsafe_execute!(plan::cFFTWPlan{T},
                X::StridedArray{T}, Y::StridedArray{T}) where {T<:fftwSingle} =
    ccall((:fftwf_execute_dft,libfftw3f), Cvoid,
          (PlanPtr,Ptr{T},Ptr{T}), plan, X, Y)

unsafe_execute!(plan::rFFTWPlan{Float64,FORWARD},
                X::StridedArray{Float64}, Y::StridedArray{Complex{Float64}}) =
    ccall((:fftw_execute_dft_r2c,libfftw3), Cvoid,
          (PlanPtr,Ptr{Float64},Ptr{Complex{Float64}}), plan, X, Y)

unsafe_execute!(plan::rFFTWPlan{Float32,FORWARD},
                X::StridedArray{Float32}, Y::StridedArray{Complex{Float32}}) =
    ccall((:fftwf_execute_dft_r2c,libfftw3f), Cvoid,
          (PlanPtr,Ptr{Float32},Ptr{Complex{Float32}}), plan, X, Y)

unsafe_execute!(plan::rFFTWPlan{Complex{Float64},BACKWARD},
                X::StridedArray{Complex{Float64}}, Y::StridedArray{Float64}) =
    ccall((:fftw_execute_dft_c2r,libfftw3), Cvoid,
          (PlanPtr,Ptr{Complex{Float64}},Ptr{Float64}), plan, X, Y)

unsafe_execute!(plan::rFFTWPlan{Complex{Float32},BACKWARD},
                X::StridedArray{Complex{Float32}}, Y::StridedArray{Float32}) =
    ccall((:fftwf_execute_dft_c2r,libfftw3f), Cvoid,
          (PlanPtr,Ptr{Complex{Float32}},Ptr{Float32}), plan, X, Y)

unsafe_execute!(plan::r2rFFTWPlan{T},
                X::StridedArray{T}, Y::StridedArray{T}) where {T<:fftwDouble} =
    ccall((:fftw_execute_r2r,libfftw3), Cvoid,
          (PlanPtr,Ptr{T},Ptr{T}), plan, X, Y)

unsafe_execute!(plan::r2rFFTWPlan{T},
                X::StridedArray{T}, Y::StridedArray{T}) where {T<:fftwSingle} =
    ccall((:fftwf_execute_r2r,libfftw3f), Cvoid,
          (PlanPtr,Ptr{T},Ptr{T}), plan, X, Y)

# NOTE ON GC (garbage collection):
#    The FFTWPlan has a finalizer so that gc will destroy the plan,
#    which is necessary for gc to work with plan_fft.  However,
#    even when we are creating a single-use FFTWPlan [e.g. for fftn(x)],
#    we intentionally do NOT call destroy_plan explicitly, and instead
#    wait for garbage collection.  The reason is that, in the common
#    case where the user calls fft(x) a second time soon afterwards,
#    if destroy_plan has not yet been called then FFTW will internally
#    re-use the table of trigonometric constants from the first plan.

# Compute dims and howmany for FFTW guru planner
_anyrepeated(::Union{Number, AbstractUnitRange}) = false
function _anyrepeated(region)
    any(region) do x
        count(==(x), region) > 1
    end
end

# Utility methods to reduce allocations in dims_howmany
@inline _setindex(oreg, v, n) = (oreg[n] = v; oreg)
@inline _setindex(oreg::Tuple, v, n) = Base.setindex(oreg, v, n)
@inline _filtercoll(region::Union{Int, Tuple}, len) = ntuple(zero, len)
@inline _filtercoll(region, len) = Vector{Int}(undef, len)
# Optimized filter(∉(region), 1:ndims(X))
function _filter_notin_region(region, ::Val{ndimsX}) where {ndimsX}
    oreg = _filtercoll(region, ndimsX - length(region))
    n = 1
    for dim in 1:ndimsX
        dim in region && continue
        oreg = _setindex(oreg, dim, n)
        n += 1
    end
    oreg
end
function dims_howmany(X::StridedArray, Y::StridedArray, sz, region)
    if _anyrepeated(region)
        throw(ArgumentError("each dimension can be transformed at most once"))
    end
    ist = strides(X)
    ost = strides(Y)
    dims = Matrix{Int}(undef, 3, length(region))
    for (ind, i) in enumerate(region)
        dims[1, ind] = sz[i]
        dims[2, ind] = ist[i]
        dims[3, ind] = ost[i]
    end

    oreg = _filter_notin_region(region, Val(ndims(X)))
    howmany = Matrix{Int}(undef, 3, length(oreg))
    for (ind, i) in enumerate(oreg)
        howmany[1, ind] = sz[i]
        howmany[2, ind] = ist[i]
        howmany[3, ind] = ost[i]
    end

    return dims, howmany
end

function fix_kinds(region, kinds)
    if length(kinds) != length(region)
        if length(kinds) > length(region)
            throw(ArgumentError("too many transform kinds"))
        else
            if isempty(kinds)
                throw(ArgumentError("must supply a transform kind"))
            end
            k = Vector{Int32}(undef, length(region))
            copyto!(k, kinds)
            k[length(kinds)+1:end] .= kinds[end]
        end
    else
        k = Vector{Int32}(undef, length(kinds))
        k .= kinds
    end
    if any(x -> x < 0 || x > 10, k)
        throw(ArgumentError("invalid transform kind"))
    end
    return k
end

_circshiftmin1(v) = circshift(collect(Int, v), -1)
_circshiftmin1(t::Tuple) = (t[2:end]..., t[1])
_circshiftmin1(x::Integer) = x

# low-level FFTWPlan creation (for internal use in FFTW module)
for (Tr,Tc,fftw,lib) in ((:Float64,:(Complex{Float64}),"fftw",:libfftw3),
                         (:Float32,:(Complex{Float32}),"fftwf",:libfftw3f))
    @eval @exclusive function cFFTWPlan{$Tc,K,inplace,N}(X::StridedArray{$Tc,N},
                                              Y::StridedArray{$Tc,N},
                                              region, flags::Integer, timelimit::Real) where {K,inplace,N}
        direction = K
        unsafe_set_timelimit($Tr, timelimit)
        R = isa(region, Tuple) ? region : copy(region)
        dims, howmany = dims_howmany(X, Y, size(X), R)
        plan = ccall(($(string(fftw,"_plan_guru64_dft")),$lib),
                     PlanPtr,
                     (Int32, Ptr{Int}, Int32, Ptr{Int},
                      Ptr{$Tc}, Ptr{$Tc}, Int32, UInt32),
                     size(dims,2), dims, size(howmany,2), howmany,
                     X, Y, direction, flags)
        unsafe_set_timelimit($Tr, NO_TIMELIMIT)
        if plan == C_NULL
            error("FFTW could not create plan") # shouldn't normally happen
        end
        return cFFTWPlan{$Tc,K,inplace,N}(plan, flags, R, X, Y)
    end

    @eval @exclusive function rFFTWPlan{$Tr,$FORWARD,inplace,N}(X::StridedArray{$Tr,N},
                                                     Y::StridedArray{$Tc,N},
                                                     region, flags::Integer, timelimit::Real) where {inplace,N}
        R = isa(region, Tuple) ? region : copy(region)
        regionshft = _circshiftmin1(region) # FFTW halves last dim
        unsafe_set_timelimit($Tr, timelimit)
        dims, howmany = dims_howmany(X, Y, size(X), regionshft)
        plan = ccall(($(string(fftw,"_plan_guru64_dft_r2c")),$lib),
                     PlanPtr,
                     (Int32, Ptr{Int}, Int32, Ptr{Int},
                      Ptr{$Tr}, Ptr{$Tc}, UInt32),
                     size(dims,2), dims, size(howmany,2), howmany,
                     X, Y, flags)
        unsafe_set_timelimit($Tr, NO_TIMELIMIT)
        if plan == C_NULL
            error("FFTW could not create plan") # shouldn't normally happen
        end
        return rFFTWPlan{$Tr,$FORWARD,inplace,N}(plan, flags, R, X, Y)
    end

    @eval @exclusive function rFFTWPlan{$Tc,$BACKWARD,inplace,N}(X::StridedArray{$Tc,N},
                                                      Y::StridedArray{$Tr,N},
                                                      region, flags::Integer, timelimit::Real) where {inplace,N}
        R = isa(region, Tuple) ? region : copy(region)
        regionshft = _circshiftmin1(region) # FFTW halves last dim
        unsafe_set_timelimit($Tr, timelimit)
        dims, howmany = dims_howmany(X, Y, size(Y), regionshft)
        plan = ccall(($(string(fftw,"_plan_guru64_dft_c2r")),$lib),
                     PlanPtr,
                     (Int32, Ptr{Int}, Int32, Ptr{Int},
                      Ptr{$Tc}, Ptr{$Tr}, UInt32),
                     size(dims,2), dims, size(howmany,2), howmany,
                     X, Y, flags)
        unsafe_set_timelimit($Tr, NO_TIMELIMIT)
        if plan == C_NULL
            error("FFTW could not create plan") # shouldn't normally happen
        end
        return rFFTWPlan{$Tc,$BACKWARD,inplace,N}(plan, flags, R, X, Y)
    end

    @eval @exclusive function r2rFFTWPlan{$Tr,Any,inplace,N}(X::StridedArray{$Tr,N},
                                                  Y::StridedArray{$Tr,N},
                                                  region, kinds, flags::Integer,
                                                  timelimit::Real) where {inplace,N}

        R = isa(region, Tuple) ? region : copy(region)
        knd = fix_kinds(region, kinds)
        unsafe_set_timelimit($Tr, timelimit)
        dims, howmany = dims_howmany(X, Y, size(X), region)
        plan = ccall(($(string(fftw,"_plan_guru64_r2r")),$lib),
                     PlanPtr,
                     (Int32, Ptr{Int}, Int32, Ptr{Int},
                      Ptr{$Tr}, Ptr{$Tr}, Ptr{Int32}, UInt32),
                     size(dims,2), dims, size(howmany,2), howmany,
                     X, Y, knd, flags)
        unsafe_set_timelimit($Tr, NO_TIMELIMIT)
        if plan == C_NULL
            error("FFTW could not create plan") # shouldn't normally happen
        end
        r2rFFTWPlan{$Tr,typeof(knd),inplace,N}(plan, flags, R, X, Y, knd)
    end

    # support r2r transforms of complex = transforms of real & imag parts
    @eval @exclusive function r2rFFTWPlan{$Tc,Any,inplace,N}(X::StridedArray{$Tc,N},
                                                  Y::StridedArray{$Tc,N},
                                                  region, kinds, flags::Integer,
                                                  timelimit::Real) where {inplace,N}

        R = isa(region, Tuple) ? region : copy(region)
        knd = fix_kinds(region, kinds)
        unsafe_set_timelimit($Tr, timelimit)
        dims, howmany = dims_howmany(X, Y, size(X), region)
        @views begin
            dims[2:3, :] .*= 2
            howmany[2:3, :] .*= 2
        end
        howmany = [howmany [2,1,1]] # append loop over real/imag parts
        plan = ccall(($(string(fftw,"_plan_guru64_r2r")),$lib),
                     PlanPtr,
                     (Int32, Ptr{Int}, Int32, Ptr{Int},
                      Ptr{$Tc}, Ptr{$Tc}, Ptr{Int32}, UInt32),
                     size(dims,2), dims, size(howmany,2), howmany,
                     X, Y, knd, flags)
        unsafe_set_timelimit($Tr, NO_TIMELIMIT)
        if plan == C_NULL
            error("FFTW could not create plan") # shouldn't normally happen
        end
        r2rFFTWPlan{$Tc,typeof(knd),inplace,N}(plan, flags, R, X, Y, knd)
    end
end

# Convert arrays of numeric types to FFTW-supported packed complex-float types
# (FIXME: is there a way to use the Julia promotion rules more cleverly here?)
fftwcomplex(X::StridedArray{<:fftwComplex}) = X
fftwcomplex(X::AbstractArray{T}) where {T<:fftwReal} =
    copyto!(Array{typeof(complex(zero(T)))}(undef, size(X)), X)
fftwcomplex(X::AbstractArray{<:Real}) = copyto!(Array{Complex{Float64}}(undef, size(X)),X)
fftwcomplex(X::AbstractArray{<:Complex}) = copyto!(Array{Complex{Float64}}(undef, size(X)), X)
fftwfloat(X::StridedArray{<:fftwReal}) = X
fftwfloat(X::AbstractArray{<:Real}) = copyto!(Array{Float64}(undef, size(X)), X)
fftwfloat(X::AbstractArray{<:Complex}) = fftwcomplex(X)

for (f,direction) in ((:fft,FORWARD), (:bfft,BACKWARD))
    plan_f = Symbol("plan_",f)
    plan_f! = Symbol("plan_",f,"!")
    idirection = -direction
    @eval begin
        function $plan_f(X::StridedArray{T,N}, region;
                         flags::Integer=ESTIMATE,
                         timelimit::Real=NO_TIMELIMIT,
                         num_threads::Union{Nothing, Integer} = nothing) where {T<:fftwComplex,N}
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    $plan_f(X, region; flags = flags, timelimit = timelimit)
                end
                return plan
            end
            cFFTWPlan{T,$direction,false,N}(X, fakesimilar(flags, X, T),
                                            region, flags, timelimit)
        end

        function $plan_f!(X::StridedArray{T,N}, region;
                          flags::Integer=ESTIMATE,
                          timelimit::Real=NO_TIMELIMIT,
                          num_threads::Union{Nothing, Integer} = nothing ) where {T<:fftwComplex,N}
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    $plan_f!(X, region; flags = flags, timelimit = timelimit)
                end
                return plan
            end
            cFFTWPlan{T,$direction,true,N}(X, X, region, flags, timelimit)
        end
        $plan_f(X::StridedArray{<:fftwComplex}; kws...) =
            $plan_f(X, ntuple(identity, ndims(X)); kws...)
        $plan_f!(X::StridedArray{<:fftwComplex}; kws...) =
            $plan_f!(X, ntuple(identity, ndims(X)); kws...)

        function plan_inv(p::cFFTWPlan{T,$direction,inplace,N};
                          num_threads::Union{Nothing, Integer} = nothing) where {T<:fftwComplex,N,inplace}
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    plan_inv(p)
                end
                return plan
            end
            X = Array{T}(undef, p.sz)
            Y = inplace ? X : fakesimilar(p.flags, X, T)
            ScaledPlan(cFFTWPlan{T,$idirection,inplace,N}(X, Y, p.region,
                                                          p.flags, NO_TIMELIMIT),
                       normalization(X, p.region))
        end
    end
end

function mul!(y::StridedArray{T}, p::cFFTWPlan{T}, x::StridedArray{T}) where T
    assert_applicable(p, x, y)
    unsafe_execute!(p, x, y)
    return y
end

function *(p::cFFTWPlan{T,K,false}, x::StridedArray{T,N}) where {T,K,N}
    assert_applicable(p, x)
    y = Array{T}(undef, p.osz)::Array{T,N}
    unsafe_execute!(p, x, y)
    return y
end

function *(p::cFFTWPlan{T,K,true}, x::StridedArray{T}) where {T,K}
    assert_applicable(p, x)
    unsafe_execute!(p, x, x)
    return x
end

# rfft/brfft and planned variants.  No in-place version for now.

for (Tr,Tc) in ((:Float32,:(Complex{Float32})),(:Float64,:(Complex{Float64})))
    # Note: use $FORWARD and $BACKWARD below because of issue #9775
    @eval begin
        function plan_rfft(X::StridedArray{$Tr,N}, region;
                           flags::Integer=ESTIMATE,
                           timelimit::Real=NO_TIMELIMIT,
                           num_threads::Union{Nothing, Integer} = nothing) where N
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    plan_rfft(X, region; flags = flags, timelimit = timelimit)
                end
                return plan
            end
            osize = rfft_output_size(X, region)
            Y = flags&ESTIMATE != 0 ? FakeArray{$Tc}(osize) : Array{$Tc}(undef, osize)
            rFFTWPlan{$Tr,$FORWARD,false,N}(X, Y, region, flags, timelimit)
        end

        function plan_brfft(X::StridedArray{$Tc,N}, d::Integer, region;
                            flags::Integer=ESTIMATE,
                            timelimit::Real=NO_TIMELIMIT,
                            num_threads::Union{Nothing, Integer} = nothing) where N
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    plan_brfft(X, d, region; flags = flags, timelimit = timelimit)
                end
                return plan
            end
            osize = brfft_output_size(X, d, region)
            Y = flags&ESTIMATE != 0 ? FakeArray{$Tr}(osize) : Array{$Tr}(undef, osize)

            # FFTW currently doesn't support PRESERVE_INPUT for
            # multidimensional out-of-place c2r transforms, so
            # we have to handle 1d and >1d cases separately with a copy.  Ugh.
            if length(region) <= 1
                rFFTWPlan{$Tc,$BACKWARD,false,N}(X, Y, region,
                                                 flags | PRESERVE_INPUT,
                                                 timelimit)
            else
                rFFTWPlan{$Tc,$BACKWARD,false,N}(copy(X), Y, region, flags,
                                                 timelimit)
            end
        end

        plan_rfft(X::StridedArray{$Tr};kws...)=plan_rfft(X,ntuple(identity, ndims(X));kws...)
        plan_brfft(X::StridedArray{$Tr};kws...)=plan_brfft(X,ntuple(identity, ndims(X));kws...)

        function plan_inv(p::rFFTWPlan{$Tr,$FORWARD,false,N},
                          num_threads::Union{Nothing, Integer} = nothing) where N
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    plan_inv(p)
                end
                return plan
            end
            X = Array{$Tr}(undef, p.sz)
            Y = p.flags&ESTIMATE != 0 ? FakeArray{$Tc}(p.osz) : Array{$Tc}(undef, p.osz)
            ScaledPlan(rFFTWPlan{$Tc,$BACKWARD,false,N}(Y, X, p.region,
                                                        length(p.region) <= 1 ?
                                                        p.flags | PRESERVE_INPUT :
                                                        p.flags, NO_TIMELIMIT),
                       normalization(X, p.region))
        end

        function plan_inv(p::rFFTWPlan{$Tc,$BACKWARD,false,N};
                          num_threads::Union{Nothing, Integer} = nothing) where N
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    plan_inv(p)
                end
                return plan
            end
            X = Array{$Tc}(undef, p.sz)
            Y = p.flags&ESTIMATE != 0 ? FakeArray{$Tr}(p.osz) : Array{$Tr}(undef, p.osz)
            ScaledPlan(rFFTWPlan{$Tr,$FORWARD,false,N}(Y, X, p.region,
                                                       p.flags, NO_TIMELIMIT),
                       normalization(Y, p.region))
        end

        function mul!(y::StridedArray{$Tc}, p::rFFTWPlan{$Tr,$FORWARD}, x::StridedArray{$Tr})
            assert_applicable(p, x, y)
            unsafe_execute!(p, x, y)
            return y
        end
        function mul!(y::StridedArray{$Tr}, p::rFFTWPlan{$Tc,$BACKWARD}, x::StridedArray{$Tc})
            assert_applicable(p, x, y)
            unsafe_execute!(p, x, y) # note: may overwrite x as well as y!
            return y
        end

        function *(p::rFFTWPlan{$Tr,$FORWARD,false}, x::StridedArray{$Tr,N}) where N
            assert_applicable(p, x)
            y = Array{$Tc}(undef, p.osz)
            unsafe_execute!(p, x, y)
            return y
        end

        function *(p::rFFTWPlan{$Tc,$BACKWARD,false}, x::StridedArray{$Tc,N}) where N
            if p.flags & PRESERVE_INPUT != 0
                assert_applicable(p, x)
                y = Array{$Tr}(undef, p.osz)
                unsafe_execute!(p, x, y)
            else # need to make a copy to avoid overwriting x
                xc = copy(x)
                assert_applicable(p, xc)
                y = Array{$Tr}(undef, p.osz)
                unsafe_execute!(p, xc, y)
            end
            return y
        end
    end
end

# FFTW r2r transforms (low-level interface)

_ntupleid(v) = ntuple(identity, v)
for f in (:r2r, :r2r!)
    pf = Symbol("plan_", f)
    @eval begin
        $f(x::AbstractArray{<:fftwNumber}, kinds) = $pf(x, kinds) * x
        $f(x::AbstractArray{<:fftwNumber}, kinds, region) = $pf(x, kinds, region) * x
        $pf(x::AbstractArray, kinds; kws...) = $pf(x, kinds, _ntupleid(Val(ndims(x))); kws...)
        $f(x::AbstractArray{<:Real}, kinds, region=_ntupleid(Val(ndims(x)))) = $f(fftwfloat(x), kinds, region)
        $pf(x::AbstractArray{<:Real}, kinds, region; kws...) = $pf(fftwfloat(x), kinds, region; kws...)
        $f(x::AbstractArray{<:Complex}, kinds, region=_ntupleid(Val(ndims(x)))) = $f(fftwcomplex(x), kinds, region)
        $pf(x::AbstractArray{<:Complex}, kinds, region; kws...) = $pf(fftwcomplex(x), kinds, region; kws...)
    end
end

function plan_r2r(X::StridedArray{T,N}, kinds, region;
                  flags::Integer=ESTIMATE,
                  timelimit::Real=NO_TIMELIMIT,
                  num_threads::Union{Nothing, Integer} = nothing) where {T<:fftwNumber,N}
    if num_threads !== nothing
        plan = set_num_threads(num_threads) do
            plan_r2r(X, kinds, region; flags = flags, timelimit = timelimit)
        end
        return plan
    end
    r2rFFTWPlan{T,Any,false,N}(X, fakesimilar(flags, X, T), region, kinds,
                               flags, timelimit)
end

function plan_r2r!(X::StridedArray{T,N}, kinds, region;
                   flags::Integer=ESTIMATE,
                   timelimit::Real=NO_TIMELIMIT,
                   num_threads::Union{Nothing, Integer} = nothing) where {T<:fftwNumber,N}
    if num_threads !== nothing
        plan = set_num_threads(num_threads) do
            plan_r2r!(X, kinds, region; flags = flags, timelimit = timelimit)
        end
        return plan
    end
    r2rFFTWPlan{T,Any,true,N}(X, X, region, kinds, flags, timelimit)
end

# mapping from r2r kind to the corresponding inverse transform
const inv_kind = Dict{Int,Int}(R2HC => HC2R, HC2R => R2HC, DHT => DHT,
                               REDFT00 => REDFT00,
                               REDFT01 => REDFT10, REDFT10 => REDFT01,
                               REDFT11 => REDFT11,
                               RODFT00 => RODFT00,
                               RODFT01 => RODFT10, RODFT10 => RODFT01,
                               RODFT11 => RODFT11)

# r2r inverses are normalized to 1/N, where N is a "logical" size
# the transform with length n and kind k:
function logical_size(n::Integer, k::Integer)
    k <= DHT && return n
    k == REDFT00 && return 2(n-1)
    k == RODFT00 && return 2(n+1)
    return 2n
end

function plan_inv(p::r2rFFTWPlan{T,K,inplace,N};
                  num_threads::Union{Nothing, Integer} = nothing) where {T<:fftwNumber,K,inplace,N}
    if num_threads !== nothing
        set_num_threads(num_threads) do
            plan = plan_inv(p)
        end
        return plan
    end
    X = Array{T}(undef, p.sz)
    # broadcast getindex to preserve tuples
    iK = fix_kinds(p.region, getindex.((inv_kind,), p.kinds))
    Y = inplace ? X : fakesimilar(p.flags, X, T)
    ScaledPlan(r2rFFTWPlan{T,Any,inplace,N}(X, Y, p.region, iK,
                                            p.flags, NO_TIMELIMIT),
               normalization(real(T),
                             map(logical_size, [p.sz...][[p.region...]], iK),
                             1:length(p.region)))
end

function mul!(y::StridedArray{T}, p::r2rFFTWPlan{T}, x::StridedArray{T}) where T
    assert_applicable(p, x, y)
    unsafe_execute!(p, x, y)
    return y
end

function *(p::r2rFFTWPlan{T,K,false}, x::StridedArray{T,N}) where {T,K,N}
    assert_applicable(p, x)
    y = Array{T}(undef, p.osz)::Array{T,N}
    unsafe_execute!(p, x, y)
    return y
end

function *(p::r2rFFTWPlan{T,K,true}, x::StridedArray{T}) where {T,K}
    assert_applicable(p, x)
    unsafe_execute!(p, x, x)
    return x
end

#######################################################################

AbstractFFTs.AdjointStyle(::cFFTWPlan) = AbstractFFTs.FFTAdjointStyle()
AbstractFFTs.AdjointStyle(::rFFTWPlan{T, FORWARD}) where {T} = AbstractFFTs.RFFTAdjointStyle()
AbstractFFTs.AdjointStyle(p::rFFTWPlan{T, BACKWARD}) where {T} = AbstractFFTs.IRFFTAdjointStyle(p.osz[first(p.region)])
