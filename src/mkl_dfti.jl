# MKL native DFTI interface for FFTW.jl
# Uses MKL's DftiCreateDescriptor/DftiComputeForward/DftiComputeBackward API
# instead of the limited FFTW3 compatibility wrappers (guru64), which fail
# for howmany_rank > 1 (e.g., transforming one dim of a 3D+ array).

# This file is only included when fftw_provider == "mkl".

#==============================================================================#
# MKL_LONG type: C `long`, which is 64-bit on Linux/macOS but 32-bit on Windows
#==============================================================================#
const MKL_LONG = Clong

#==============================================================================#
# DFTI configuration constants (from mkl_dfti.h)
#==============================================================================#

# DFTI_CONFIG_VALUE enum
const DFTI_COMMITTED       = MKL_LONG(30)
const DFTI_UNCOMMITTED     = MKL_LONG(31)

# Precision
const DFTI_SINGLE          = MKL_LONG(35)
const DFTI_DOUBLE          = MKL_LONG(36)

# Forward domain
const DFTI_COMPLEX         = MKL_LONG(32)
const DFTI_REAL            = MKL_LONG(33)

# Placement
const DFTI_INPLACE         = MKL_LONG(43)
const DFTI_NOT_INPLACE     = MKL_LONG(44)

# Complex storage
const DFTI_COMPLEX_COMPLEX = MKL_LONG(39)
const DFTI_REAL_REAL       = MKL_LONG(40)
const DFTI_COMPLEX_REAL    = MKL_LONG(41)

# Packed format
const DFTI_CCE_FORMAT      = MKL_LONG(56)
const DFTI_CCS_FORMAT      = MKL_LONG(45)
const DFTI_PACK_FORMAT     = MKL_LONG(46)
const DFTI_PERM_FORMAT     = MKL_LONG(47)

# Allow/Avoid
const DFTI_ALLOW           = MKL_LONG(51)
const DFTI_AVOID           = MKL_LONG(52)

# DFTI_CONFIG_PARAM enum
const DFTI_FORWARD_DOMAIN           = MKL_LONG(0)
const DFTI_DIMENSION                = MKL_LONG(1)
const DFTI_LENGTHS                  = MKL_LONG(2)
const DFTI_PRECISION                = MKL_LONG(3)
const DFTI_FORWARD_SCALE            = MKL_LONG(4)
const DFTI_BACKWARD_SCALE           = MKL_LONG(5)
const DFTI_NUMBER_OF_TRANSFORMS     = MKL_LONG(7)
const DFTI_COMPLEX_STORAGE          = MKL_LONG(8)
const DFTI_REAL_STORAGE             = MKL_LONG(9)
const DFTI_CONJUGATE_EVEN_STORAGE   = MKL_LONG(10)
const DFTI_PLACEMENT                = MKL_LONG(11)
const DFTI_INPUT_STRIDES            = MKL_LONG(12)
const DFTI_OUTPUT_STRIDES           = MKL_LONG(13)
const DFTI_INPUT_DISTANCE           = MKL_LONG(14)
const DFTI_OUTPUT_DISTANCE          = MKL_LONG(15)
const DFTI_WORKSPACE                = MKL_LONG(17)
const DFTI_ORDERING                 = MKL_LONG(18)
const DFTI_TRANSPOSE                = MKL_LONG(19)
const DFTI_DESCRIPTOR_NAME          = MKL_LONG(20)
const DFTI_PACKED_FORMAT            = MKL_LONG(21)
const DFTI_COMMIT_STATUS            = MKL_LONG(22)
const DFTI_VERSION                  = MKL_LONG(23)
const DFTI_NUMBER_OF_USER_THREADS   = MKL_LONG(26)
const DFTI_THREAD_LIMIT             = MKL_LONG(27)
const DFTI_DESTROY_INPUT            = MKL_LONG(28)

# Ordered/Scrambled
const DFTI_ORDERED             = MKL_LONG(48)
const DFTI_BACKWARD_SCRAMBLED  = MKL_LONG(49)

#==============================================================================#
# DFTI API wrappers
#==============================================================================#

const DftiDescriptor = Ptr{Cvoid}

# Note: DftiCreateDescriptor is a macro in mkl_dfti.h that dispatches to
# specific functions based on precision and dimensionality. We call the
# specific underlying functions directly.

function dfti_create_descriptor_1d(precision::MKL_LONG, domain::MKL_LONG, length::MKL_LONG)
    handle = Ref{DftiDescriptor}(C_NULL)
    if precision == DFTI_SINGLE
        status = ccall((:DftiCreateDescriptor_s_1d, libfftw3), MKL_LONG,
                       (Ref{DftiDescriptor}, MKL_LONG, MKL_LONG),
                       handle, domain, length)
    else
        status = ccall((:DftiCreateDescriptor_d_1d, libfftw3), MKL_LONG,
                       (Ref{DftiDescriptor}, MKL_LONG, MKL_LONG),
                       handle, domain, length)
    end
    status != 0 && error("DftiCreateDescriptor failed: $(dfti_error_message(status))")
    return handle[]
end

function dfti_create_descriptor_md(precision::MKL_LONG, domain::MKL_LONG, ndim::MKL_LONG, lengths::Vector{MKL_LONG})
    handle = Ref{DftiDescriptor}(C_NULL)
    if precision == DFTI_SINGLE
        status = ccall((:DftiCreateDescriptor_s_md, libfftw3), MKL_LONG,
                       (Ref{DftiDescriptor}, MKL_LONG, MKL_LONG, Ptr{MKL_LONG}),
                       handle, domain, ndim, lengths)
    else
        status = ccall((:DftiCreateDescriptor_d_md, libfftw3), MKL_LONG,
                       (Ref{DftiDescriptor}, MKL_LONG, MKL_LONG, Ptr{MKL_LONG}),
                       handle, domain, ndim, lengths)
    end
    status != 0 && error("DftiCreateDescriptor failed: $(dfti_error_message(status))")
    return handle[]
end

function dfti_create_descriptor(precision::MKL_LONG, domain::MKL_LONG, lengths::Vector{MKL_LONG})
    ndim = length(lengths)
    if ndim == 1
        return dfti_create_descriptor_1d(precision, domain, lengths[1])
    else
        return dfti_create_descriptor_md(precision, domain, MKL_LONG(ndim), lengths)
    end
end

function dfti_set_value(handle::DftiDescriptor, param::MKL_LONG, value::MKL_LONG)
    status = ccall((:DftiSetValue, libfftw3), MKL_LONG,
                   (DftiDescriptor, MKL_LONG, MKL_LONG),
                   handle, param, value)
    status != 0 && error("DftiSetValue failed (param=$param): $(dfti_error_message(status))")
end

function dfti_set_value(handle::DftiDescriptor, param::MKL_LONG, value::Float32)
    status = ccall((:DftiSetValue, libfftw3), MKL_LONG,
                   (DftiDescriptor, MKL_LONG, Float32),
                   handle, param, value)
    status != 0 && error("DftiSetValue failed (param=$param): $(dfti_error_message(status))")
end

function dfti_set_value(handle::DftiDescriptor, param::MKL_LONG, value::Float64)
    status = ccall((:DftiSetValue, libfftw3), MKL_LONG,
                   (DftiDescriptor, MKL_LONG, Float64),
                   handle, param, value)
    status != 0 && error("DftiSetValue failed (param=$param): $(dfti_error_message(status))")
end

function dfti_set_value(handle::DftiDescriptor, param::MKL_LONG, value::Vector{MKL_LONG})
    status = ccall((:DftiSetValue, libfftw3), MKL_LONG,
                   (DftiDescriptor, MKL_LONG, Ptr{MKL_LONG}),
                   handle, param, value)
    status != 0 && error("DftiSetValue failed (param=$param): $(dfti_error_message(status))")
end

function dfti_commit_descriptor(handle::DftiDescriptor)
    status = ccall((:DftiCommitDescriptor, libfftw3), MKL_LONG,
                   (DftiDescriptor,), handle)
    status != 0 && error("DftiCommitDescriptor failed: $(dfti_error_message(status))")
end

function dfti_compute_forward_inplace(handle::DftiDescriptor, data::Ptr)
    status = ccall((:DftiComputeForward, libfftw3), MKL_LONG,
                   (DftiDescriptor, Ptr{Cvoid}),
                   handle, data)
    status != 0 && error("DftiComputeForward failed: $(dfti_error_message(status))")
end

function dfti_compute_forward_outofplace(handle::DftiDescriptor, indata::Ptr, outdata::Ptr)
    status = ccall((:DftiComputeForward, libfftw3), MKL_LONG,
                   (DftiDescriptor, Ptr{Cvoid}, Ptr{Cvoid}),
                   handle, indata, outdata)
    status != 0 && error("DftiComputeForward failed: $(dfti_error_message(status))")
end

function dfti_compute_backward_inplace(handle::DftiDescriptor, data::Ptr)
    status = ccall((:DftiComputeBackward, libfftw3), MKL_LONG,
                   (DftiDescriptor, Ptr{Cvoid}),
                   handle, data)
    status != 0 && error("DftiComputeBackward failed: $(dfti_error_message(status))")
end

function dfti_compute_backward_outofplace(handle::DftiDescriptor, indata::Ptr, outdata::Ptr)
    status = ccall((:DftiComputeBackward, libfftw3), MKL_LONG,
                   (DftiDescriptor, Ptr{Cvoid}, Ptr{Cvoid}),
                   handle, indata, outdata)
    status != 0 && error("DftiComputeBackward failed: $(dfti_error_message(status))")
end

function dfti_free_descriptor(handle::DftiDescriptor)
    ref = Ref(handle)
    ccall((:DftiFreeDescriptor, libfftw3), MKL_LONG,
          (Ref{DftiDescriptor},), ref)
end

function dfti_error_message(status::MKL_LONG)
    ptr = ccall((:DftiErrorMessage, libfftw3), Ptr{UInt8}, (MKL_LONG,), status)
    ptr == C_NULL ? "unknown error (status=$status)" : unsafe_string(ptr)
end

#==============================================================================#
# MKL Plan types
#==============================================================================#

# MKL DFTI plan for complex-to-complex transforms
mutable struct MKLcPlan{T<:fftwComplex,K,inplace,N,G} <: FFTWPlan{T,K,inplace}
    handle::DftiDescriptor
    sz::NTuple{N,Int}
    osz::NTuple{N,Int}
    istride::NTuple{N,Int}
    ostride::NTuple{N,Int}
    ialign::Int32
    oalign::Int32
    flags::UInt32
    region::G
    outer_offsets::Vector{Tuple{Int,Int}} # (ioffset, ooffset) for outer batch dims
    pinv::ScaledPlan
    function MKLcPlan{T,K,inplace,N,G}(handle::DftiDescriptor, flags::Integer, R::G,
                                       X::StridedArray{T,N}, Y::StridedArray,
                                       outer_offsets::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[]) where {T<:fftwComplex,K,inplace,N,G}
        p = new(handle, size(X), size(Y), strides(X), strides(Y),
                alignment_of(X), alignment_of(Y), UInt32(flags), R, outer_offsets)
        finalizer(maybe_destroy_plan, p)
        p
    end
end

function MKLcPlan{T,K,inplace,N}(handle::DftiDescriptor, flags::Integer, R::G,
                                 X::StridedArray{T,N}, Y::StridedArray,
                                 outer_offsets::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[]) where {T<:fftwComplex,K,inplace,N,G}
    MKLcPlan{T,K,inplace,N,G}(handle, flags, R, X, Y, outer_offsets)
end

# MKL DFTI plan for real-to-complex (FORWARD) / complex-to-real (BACKWARD) transforms
mutable struct MKLrPlan{T<:fftwNumber,K,inplace,N,G} <: FFTWPlan{T,K,inplace}
    handle::DftiDescriptor
    sz::NTuple{N,Int}
    osz::NTuple{N,Int}
    istride::NTuple{N,Int}
    ostride::NTuple{N,Int}
    ialign::Int32
    oalign::Int32
    flags::UInt32
    region::G
    outer_offsets::Vector{Tuple{Int,Int}} # (ioffset, ooffset) for outer batch dims
    pinv::ScaledPlan
    function MKLrPlan{T,K,inplace,N,G}(handle::DftiDescriptor, flags::Integer, R::G,
                                       X::StridedArray{T,N}, Y::StridedArray,
                                       outer_offsets::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[]) where {T<:fftwNumber,K,inplace,N,G}
        p = new(handle, size(X), size(Y), strides(X), strides(Y),
                alignment_of(X), alignment_of(Y), UInt32(flags), R, outer_offsets)
        finalizer(maybe_destroy_plan, p)
        p
    end
end

function MKLrPlan{T,K,inplace,N}(handle::DftiDescriptor, flags::Integer, R::G,
                                 X::StridedArray{T,N}, Y::StridedArray,
                                 outer_offsets::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[]) where {T<:fftwNumber,K,inplace,N,G}
    MKLrPlan{T,K,inplace,N,G}(handle, flags, R, X, Y, outer_offsets)
end

#==============================================================================#
# Plan Destruction (override FFTWPlan defaults which call fftw_destroy_plan)
#==============================================================================#

unsafe_destroy_plan(@nospecialize(plan::MKLcPlan)) = begin
    if plan.handle != C_NULL
        dfti_free_descriptor(plan.handle)
        plan.handle = C_NULL
    end
end

unsafe_destroy_plan(@nospecialize(plan::MKLrPlan)) = begin
    if plan.handle != C_NULL
        dfti_free_descriptor(plan.handle)
        plan.handle = C_NULL
    end
end

# Override unsafe_convert to avoid the default FFTWPlan method, which would
# access the nonexistent `p.plan` field. Instead, wrap the DFTI handle as a PlanPtr.
unsafe_convert(::Type{PlanPtr}, p::MKLcPlan) = reinterpret(PlanPtr, p.handle)
unsafe_convert(::Type{PlanPtr}, p::MKLrPlan) = reinterpret(PlanPtr, p.handle)

#==============================================================================#
# Helpers for DFTI configuration
#==============================================================================#

_dfti_precision(::Type{Float64}) = DFTI_DOUBLE
_dfti_precision(::Type{Float32}) = DFTI_SINGLE
_dfti_precision(::Type{Complex{Float64}}) = DFTI_DOUBLE
_dfti_precision(::Type{Complex{Float32}}) = DFTI_SINGLE

"""
    _collect_region(region) -> Vector{Int}

Convert region from various types (Int, Tuple, UnitRange, Vector) to a Vector{Int}.
"""
_collect_region(r::Integer) = [Int(r)]
_collect_region(r::Tuple) = collect(Int, r)
_collect_region(r::AbstractVector) = collect(Int, r)
_collect_region(r::AbstractUnitRange) = collect(Int, r)

"""
    _make_dfti_complex_plan(T, K, inplace, X, Y, region, flags) -> MKLcPlan

Create a MKL DFTI plan for complex-to-complex transform.

The DFTI API supports multi-dimensional transforms natively and handles batch dimensions
via NUMBER_OF_TRANSFORMS + INPUT_DISTANCE/OUTPUT_DISTANCE for one batch dim.
For arrays with multiple non-transform dims, we use DFTI strides to handle the full layout.
"""
function _make_dfti_complex_plan(::Type{T}, ::Val{K}, ::Val{inplace},
                                 X::StridedArray{T,N}, Y::StridedArray{T,N},
                                 region, flags::Integer) where {T<:fftwComplex,K,inplace,N}
    reg = _collect_region(region)
    R = isa(region, Tuple) ? region : copy(region)

    sz = size(X)
    ist = strides(X)
    ost = strides(Y)

    # DFTI lengths: sizes along transformed dimensions, in the order given by region
    # MKL expects lengths in row-major (C) order for multi-dim, but since we
    # specify explicit strides, the lengths just correspond to the dimension sizes.
    dfti_lengths = MKL_LONG[sz[d] for d in reg]

    prec = _dfti_precision(T)

    handle = dfti_create_descriptor(prec, DFTI_COMPLEX, dfti_lengths)

    # Set placement
    if inplace
        dfti_set_value(handle, DFTI_PLACEMENT, DFTI_INPLACE)
    else
        dfti_set_value(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    end

    # Compute strides for DFTI
    # INPUT_STRIDES/OUTPUT_STRIDES always refer to the actual function call's
    # input/output, not fixed domains. Since we create separate descriptors
    # for forward and backward, and X is always the first arg (input) and
    # Y is always the second arg (output) to DftiCompute*, we always set:
    #   INPUT_STRIDES  = X strides (first arg to DftiCompute*)
    #   OUTPUT_STRIDES = Y strides (second arg to DftiCompute*)
    dfti_xstrides = MKL_LONG[0; [ist[d] for d in reg]]
    dfti_ystrides = MKL_LONG[0; [ost[d] for d in reg]]
    dfti_set_value(handle, DFTI_INPUT_STRIDES, dfti_xstrides)
    dfti_set_value(handle, DFTI_OUTPUT_STRIDES, dfti_ystrides)

    # Handle batch dimensions (non-transform dimensions)
    # Find all dimensions not in region
    batch_dims = [d for d in 1:N if !(d in reg)]
    outer_offsets = Tuple{Int,Int}[]

    if !isempty(batch_dims)
        # We can use NUMBER_OF_TRANSFORMS with DISTANCE for batched execution.
        # For multiple batch dims, compute total batch count and use the
        # appropriate distance. The key insight is that for column-major arrays
        # with contiguous batch dims, we can set up a single batched call.
        # Sort batch dims by input stride
        sorted_batch = sort(batch_dims, by=d->ist[d])

        min_ist_dim = sorted_batch[1]
        idist = ist[min_ist_dim]
        odist = ost[min_ist_dim]

        # Find the contiguous prefix of batch dims that can be described
        # by a single NUMBER_OF_TRANSFORMS + DISTANCE.
        expected_ist = idist
        expected_ost = odist
        inner_end = 0
        for (i, d) in enumerate(sorted_batch)
            if ist[d] != expected_ist || ost[d] != expected_ost
                break
            end
            inner_end = i
            expected_ist *= sz[d]
            expected_ost *= sz[d]
        end

        inner_count = prod(sz[sorted_batch[i]] for i in 1:inner_end)
        dfti_set_value(handle, DFTI_NUMBER_OF_TRANSFORMS, MKL_LONG(inner_count))
        dfti_set_value(handle, DFTI_INPUT_DISTANCE, MKL_LONG(idist))
        dfti_set_value(handle, DFTI_OUTPUT_DISTANCE, MKL_LONG(odist))

        # Compute outer batch offsets for dims not covered by DISTANCE
        outer_dims = sorted_batch[inner_end+1:end]
        if !isempty(outer_dims)
            offsets = Tuple{Int,Int}[(0, 0)]
            for d in outer_dims
                new_offsets = Tuple{Int,Int}[]
                for idx in 0:sz[d]-1
                    for (io, oo) in offsets
                        push!(new_offsets, (io + idx * ist[d], oo + idx * ost[d]))
                    end
                end
                offsets = new_offsets
            end
            outer_offsets = offsets
        end
    end

    dfti_commit_descriptor(handle)

    return MKLcPlan{T,K,inplace,N}(handle, flags, R, X, Y, outer_offsets)
end

"""
    _make_dfti_r2c_plan(Tr, Tc, inplace, X, Y, region, flags) -> MKLrPlan

Create a MKL DFTI plan for real-to-complex (forward rfft) transform.
"""
function _make_dfti_r2c_plan(::Type{Tr}, ::Type{Tc}, ::Val{inplace},
                              X::StridedArray{Tr,N}, Y::StridedArray{Tc,N},
                              region, flags::Integer) where {Tr<:fftwReal,Tc<:fftwComplex,inplace,N}
    reg = _collect_region(region)
    R = isa(region, Tuple) ? region : copy(region)

    sz = size(X)
    osz = size(Y)
    ist = strides(X)
    ost = strides(Y)

    # DFTI lengths are the REAL sizes along transform dimensions.
    # MKL DFTI uses C/row-major convention: the LAST element in the lengths
    # array is the innermost (fastest-varying) dimension and is the one that
    # gets halved for real transforms. Julia is column-major, so we reverse
    # the dimension order so that reg[1] (the halved dim) comes last.
    rev_reg = reverse(reg)
    dfti_lengths = MKL_LONG[sz[d] for d in rev_reg]

    prec = _dfti_precision(Tr)

    handle = dfti_create_descriptor(prec, DFTI_REAL, dfti_lengths)

    # Always out-of-place for rfft (no in-place version supported currently)
    dfti_set_value(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)

    # Use CCE (Conjugate-Complex-Even) format: output is complex array
    dfti_set_value(handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)

    # Input strides (real domain) - strides in real elements (reversed order)
    dfti_istrides = MKL_LONG[0; [ist[d] for d in rev_reg]]
    dfti_set_value(handle, DFTI_INPUT_STRIDES, dfti_istrides)

    # Output strides (complex domain) - strides in complex elements (reversed order)
    dfti_ostrides = MKL_LONG[0; [ost[d] for d in rev_reg]]
    dfti_set_value(handle, DFTI_OUTPUT_STRIDES, dfti_ostrides)

    # Handle batch dimensions
    batch_dims = [d for d in 1:N if !(d in reg)]
    outer_offsets = Tuple{Int,Int}[]
    if !isempty(batch_dims)
        sorted_batch = sort(batch_dims, by=d->ist[d])

        min_ist_dim = sorted_batch[1]
        idist = ist[min_ist_dim]
        odist = ost[min_ist_dim]

        # Find contiguous prefix of batch dims
        expected_ist = idist
        expected_ost = odist
        inner_end = 0
        for (i, d) in enumerate(sorted_batch)
            if ist[d] != expected_ist || ost[d] != expected_ost
                break
            end
            inner_end = i
            expected_ist *= sz[d]
            expected_ost *= sz[d]
        end

        inner_count = prod(sz[sorted_batch[i]] for i in 1:inner_end)
        dfti_set_value(handle, DFTI_NUMBER_OF_TRANSFORMS, MKL_LONG(inner_count))
        dfti_set_value(handle, DFTI_INPUT_DISTANCE, MKL_LONG(idist))
        dfti_set_value(handle, DFTI_OUTPUT_DISTANCE, MKL_LONG(odist))

        # Compute outer batch offsets for dims not covered by DISTANCE
        outer_dims = sorted_batch[inner_end+1:end]
        if !isempty(outer_dims)
            offsets = Tuple{Int,Int}[(0, 0)]
            for d in outer_dims
                new_offsets = Tuple{Int,Int}[]
                for idx in 0:sz[d]-1
                    for (io, oo) in offsets
                        push!(new_offsets, (io + idx * ist[d], oo + idx * ost[d]))
                    end
                end
                offsets = new_offsets
            end
            outer_offsets = offsets
        end
    end

    dfti_commit_descriptor(handle)

    return MKLrPlan{Tr,FORWARD,inplace,N}(handle, flags, R, X, Y, outer_offsets)
end

"""
    _make_dfti_c2r_plan(Tc, Tr, inplace, X, Y, region, flags) -> MKLrPlan

Create a MKL DFTI plan for complex-to-real (backward brfft) transform.
"""
function _make_dfti_c2r_plan(::Type{Tc}, ::Type{Tr}, ::Val{inplace},
                              X::StridedArray{Tc,N}, Y::StridedArray{Tr,N},
                              region, flags::Integer) where {Tc<:fftwComplex,Tr<:fftwReal,inplace,N}
    reg = _collect_region(region)
    R = isa(region, Tuple) ? region : copy(region)

    sz_in = size(X)
    sz_out = size(Y)
    ist = strides(X)
    ost = strides(Y)

    # DFTI lengths are the REAL (output) sizes along transform dimensions.
    # MKL DFTI uses C/row-major convention: reverse dims so that reg[1]
    # (the halved dim) comes last in the lengths array.
    rev_reg = reverse(reg)
    dfti_lengths = MKL_LONG[sz_out[d] for d in rev_reg]

    prec = _dfti_precision(Tr)

    handle = dfti_create_descriptor(prec, DFTI_REAL, dfti_lengths)

    # Out-of-place
    dfti_set_value(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)

    # CCE format
    dfti_set_value(handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)

    # For backward (c2r) transform, INPUT/OUTPUT refer to the backward
    # function's input/output (not the forward domain):
    #   INPUT_STRIDES  = backward input strides  = complex (X) strides
    #   OUTPUT_STRIDES = backward output strides  = real (Y) strides
    # Strides are reversed to match the reversed dimension order.
    dfti_complex_strides = MKL_LONG[0; [ist[d] for d in rev_reg]]
    dfti_real_strides = MKL_LONG[0; [ost[d] for d in rev_reg]]
    dfti_set_value(handle, DFTI_INPUT_STRIDES, dfti_complex_strides)
    dfti_set_value(handle, DFTI_OUTPUT_STRIDES, dfti_real_strides)

    # DESTROY_INPUT: MKL may overwrite input during c2r
    dfti_set_value(handle, DFTI_DESTROY_INPUT, DFTI_ALLOW)

    # Handle batch dimensions
    batch_dims = [d for d in 1:N if !(d in reg)]
    outer_offsets = Tuple{Int,Int}[]
    if !isempty(batch_dims)
        # For c2r, batch dim sizes come from the complex input array
        sorted_batch = sort(batch_dims, by=d->ist[d])

        min_ist_dim = sorted_batch[1]
        idist = ist[min_ist_dim]
        odist = ost[min_ist_dim]

        # Find contiguous prefix of batch dims
        expected_ist = idist
        expected_ost = odist
        inner_end = 0
        for (i, d) in enumerate(sorted_batch)
            if ist[d] != expected_ist || ost[d] != expected_ost
                break
            end
            inner_end = i
            expected_ist *= sz_in[d]
            expected_ost *= sz_in[d]
        end

        inner_count = prod(sz_in[sorted_batch[i]] for i in 1:inner_end)
        dfti_set_value(handle, DFTI_NUMBER_OF_TRANSFORMS, MKL_LONG(inner_count))
        dfti_set_value(handle, DFTI_INPUT_DISTANCE, MKL_LONG(idist))
        dfti_set_value(handle, DFTI_OUTPUT_DISTANCE, MKL_LONG(odist))

        # Compute outer batch offsets for dims not covered by DISTANCE
        outer_dims = sorted_batch[inner_end+1:end]
        if !isempty(outer_dims)
            offsets = Tuple{Int,Int}[(0, 0)]
            for d in outer_dims
                new_offsets = Tuple{Int,Int}[]
                for idx in 0:sz_in[d]-1
                    for (io, oo) in offsets
                        push!(new_offsets, (io + idx * ist[d], oo + idx * ost[d]))
                    end
                end
                offsets = new_offsets
            end
            outer_offsets = offsets
        end
    end

    dfti_commit_descriptor(handle)

    return MKLrPlan{Tc,BACKWARD,inplace,N}(handle, flags, R, X, Y, outer_offsets)
end

#==============================================================================#
# Execution
#==============================================================================#

function unsafe_execute!(plan::MKLcPlan{T,K,true}, X::StridedArray{T}, Y::StridedArray{T}) where {T,K}
    if isempty(plan.outer_offsets)
        if K == FORWARD
            dfti_compute_forward_inplace(plan.handle, pointer(X))
        else
            dfti_compute_backward_inplace(plan.handle, pointer(X))
        end
    else
        for (ioff, _) in plan.outer_offsets
            xp = pointer(X, ioff + 1)
            if K == FORWARD
                dfti_compute_forward_inplace(plan.handle, xp)
            else
                dfti_compute_backward_inplace(plan.handle, xp)
            end
        end
    end
end

function unsafe_execute!(plan::MKLcPlan{T,K,false}, X::StridedArray{T}, Y::StridedArray{T}) where {T,K}
    if isempty(plan.outer_offsets)
        if K == FORWARD
            dfti_compute_forward_outofplace(plan.handle, pointer(X), pointer(Y))
        else
            dfti_compute_backward_outofplace(plan.handle, pointer(X), pointer(Y))
        end
    else
        for (ioff, ooff) in plan.outer_offsets
            xp = pointer(X, ioff + 1)
            yp = pointer(Y, ooff + 1)
            if K == FORWARD
                dfti_compute_forward_outofplace(plan.handle, xp, yp)
            else
                dfti_compute_backward_outofplace(plan.handle, xp, yp)
            end
        end
    end
end

# MKLcPlan in-place execute without args (for the no-arg unsafe_execute! path)
function unsafe_execute!(plan::MKLcPlan{T,K,true}) where {T,K}
    error("MKL DFTI plans require explicit array arguments for execution")
end
function unsafe_execute!(plan::MKLcPlan{T,K,false}) where {T,K}
    error("MKL DFTI plans require explicit array arguments for execution")
end

# r2c (forward rfft)
function unsafe_execute!(plan::MKLrPlan{Tr,FORWARD}, X::StridedArray{Tr}, Y::StridedArray{<:fftwComplex}) where {Tr<:fftwReal}
    if isempty(plan.outer_offsets)
        dfti_compute_forward_outofplace(plan.handle, pointer(X), pointer(Y))
    else
        for (ioff, ooff) in plan.outer_offsets
            dfti_compute_forward_outofplace(plan.handle, pointer(X, ioff + 1), pointer(Y, ooff + 1))
        end
    end
end

# c2r (backward brfft)
function unsafe_execute!(plan::MKLrPlan{Tc,BACKWARD}, X::StridedArray{Tc}, Y::StridedArray{<:fftwReal}) where {Tc<:fftwComplex}
    if isempty(plan.outer_offsets)
        dfti_compute_backward_outofplace(plan.handle, pointer(X), pointer(Y))
    else
        for (ioff, ooff) in plan.outer_offsets
            dfti_compute_backward_outofplace(plan.handle, pointer(X, ioff + 1), pointer(Y, ooff + 1))
        end
    end
end

# No-arg versions for MKLrPlan
function unsafe_execute!(plan::MKLrPlan)
    error("MKL DFTI plans require explicit array arguments for execution")
end

#==============================================================================#
# mul! and * operators for MKL plans
#==============================================================================#

function mul!(y::StridedArray{T}, p::MKLcPlan{T}, x::StridedArray{T}) where T
    assert_applicable(p, x, y)
    unsafe_execute!(p, x, y)
    return y
end

function *(p::MKLcPlan{T,K,false}, x::StridedArray{T,N}) where {T,K,N}
    assert_applicable(p, x)
    y = Array{T}(undef, p.osz)::Array{T,N}
    unsafe_execute!(p, x, y)
    return y
end

function *(p::MKLcPlan{T,K,true}, x::StridedArray{T}) where {T,K}
    assert_applicable(p, x)
    unsafe_execute!(p, x, x)
    return x
end

for (Tr, Tc) in ((:Float32, :(Complex{Float32})), (:Float64, :(Complex{Float64})))
    @eval begin
        function mul!(y::StridedArray{$Tc}, p::MKLrPlan{$Tr,$FORWARD}, x::StridedArray{$Tr})
            assert_applicable(p, x, y)
            unsafe_execute!(p, x, y)
            return y
        end
        function mul!(y::StridedArray{$Tr}, p::MKLrPlan{$Tc,$BACKWARD}, x::StridedArray{$Tc})
            assert_applicable(p, x, y)
            unsafe_execute!(p, x, y)
            return y
        end

        function *(p::MKLrPlan{$Tr,$FORWARD,false}, x::StridedArray{$Tr,N}) where N
            assert_applicable(p, x)
            y = Array{$Tc}(undef, p.osz)
            unsafe_execute!(p, x, y)
            return y
        end

        function *(p::MKLrPlan{$Tc,$BACKWARD,false}, x::StridedArray{$Tc,N}) where N
            # c2r may overwrite input, so make a copy
            xc = copy(x)
            assert_applicable(p, xc)
            y = Array{$Tr}(undef, p.osz)
            unsafe_execute!(p, xc, y)
            return y
        end
    end
end

#==============================================================================#
# Pretty-printing
#==============================================================================#

function show(io::IO, p::MKLcPlan{T,K,inplace}) where {T,K,inplace}
    print(io, inplace ? "MKL DFTI in-place " : "MKL DFTI ",
          K < 0 ? "forward" : "backward", " plan for ")
    showfftdims(io, p.sz, p.istride, T)
end

function show(io::IO, p::MKLrPlan{T,K,inplace}) where {T,K,inplace}
    print(io, inplace ? "MKL DFTI in-place " : "MKL DFTI ",
          K < 0 ? "real-to-complex" : "complex-to-real",
          " plan for ")
    showfftdims(io, p.sz, p.istride, T)
end

#==============================================================================#
# Override plan_fft / plan_bfft for MKL
#==============================================================================#

for (f, direction) in ((:fft, FORWARD), (:bfft, BACKWARD))
    plan_f = Symbol("plan_", f)
    plan_f! = Symbol("plan_", f, "!")
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
            Y = fakesimilar(flags, X, T)
            @exclusive _make_dfti_complex_plan(T, Val($direction), Val(false), X, Y, region, flags)
        end

        function $plan_f!(X::StridedArray{T,N}, region;
                          flags::Integer=ESTIMATE,
                          timelimit::Real=NO_TIMELIMIT,
                          num_threads::Union{Nothing, Integer} = nothing) where {T<:fftwComplex,N}
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    $plan_f!(X, region; flags = flags, timelimit = timelimit)
                end
                return plan
            end
            @exclusive _make_dfti_complex_plan(T, Val($direction), Val(true), X, X, region, flags)
        end

        $plan_f(X::StridedArray{<:fftwComplex}; kws...) =
            $plan_f(X, ntuple(identity, ndims(X)); kws...)
        $plan_f!(X::StridedArray{<:fftwComplex}; kws...) =
            $plan_f!(X, ntuple(identity, ndims(X)); kws...)

        function plan_inv(p::MKLcPlan{T,$direction,inplace,N};
                          num_threads::Union{Nothing, Integer} = nothing) where {T<:fftwComplex,N,inplace}
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    plan_inv(p)
                end
                return plan
            end
            X = Array{T}(undef, p.sz)
            Y = inplace ? X : fakesimilar(p.flags, X, T)
            inv_plan = @exclusive _make_dfti_complex_plan(T, Val($idirection), Val(inplace), X, Y, p.region, p.flags)
            ScaledPlan(inv_plan, normalization(X, p.region))
        end
    end
end

#==============================================================================#
# Override plan_rfft / plan_brfft for MKL
#==============================================================================#

for (Tr, Tc) in ((:Float32, :(Complex{Float32})), (:Float64, :(Complex{Float64})))
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
            Y = flags & ESTIMATE != 0 ? FakeArray{$Tc}(osize) : Array{$Tc}(undef, osize)
            @exclusive _make_dfti_r2c_plan($Tr, $Tc, Val(false), X, Y, region, flags)
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
            Y = flags & ESTIMATE != 0 ? FakeArray{$Tr}(osize) : Array{$Tr}(undef, osize)
            @exclusive _make_dfti_c2r_plan($Tc, $Tr, Val(false), copy(X), Y, region, flags)
        end

        plan_rfft(X::StridedArray{$Tr}; kws...) = plan_rfft(X, ntuple(identity, ndims(X)); kws...)
        plan_brfft(X::StridedArray{$Tc}; kws...) = plan_brfft(X, 2*size(X,1)-1, ntuple(identity, ndims(X)); kws...)

        function plan_inv(p::MKLrPlan{$Tr,$FORWARD,false,N};
                          num_threads::Union{Nothing, Integer} = nothing) where N
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    plan_inv(p)
                end
                return plan
            end
            X = Array{$Tr}(undef, p.sz)
            Y = p.flags & ESTIMATE != 0 ? FakeArray{$Tc}(p.osz) : Array{$Tc}(undef, p.osz)
            inv_plan = @exclusive _make_dfti_c2r_plan($Tc, $Tr, Val(false), Y, X, p.region, p.flags)
            ScaledPlan(inv_plan, normalization(X, p.region))
        end

        function plan_inv(p::MKLrPlan{$Tc,$BACKWARD,false,N};
                          num_threads::Union{Nothing, Integer} = nothing) where N
            if num_threads !== nothing
                plan = set_num_threads(num_threads) do
                    plan_inv(p)
                end
                return plan
            end
            X = Array{$Tc}(undef, p.sz)
            Y = p.flags & ESTIMATE != 0 ? FakeArray{$Tr}(p.osz) : Array{$Tr}(undef, p.osz)
            inv_plan = @exclusive _make_dfti_r2c_plan($Tr, $Tc, Val(false), Y, X, p.region, p.flags)
            ScaledPlan(inv_plan, normalization(Y, p.region))
        end
    end
end

#==============================================================================#
# Adjoint styles
#==============================================================================#

AbstractFFTs.AdjointStyle(::MKLcPlan) = AbstractFFTs.FFTAdjointStyle()
AbstractFFTs.AdjointStyle(::MKLrPlan{T, FORWARD}) where {T} = AbstractFFTs.RFFTAdjointStyle()
AbstractFFTs.AdjointStyle(p::MKLrPlan{T, BACKWARD}) where {T} = AbstractFFTs.IRFFTAdjointStyle(p.osz[first(p.region)])
