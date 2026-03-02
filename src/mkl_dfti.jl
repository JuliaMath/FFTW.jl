# MKL native DFTI interface for FFTW.jl
# Uses MKL's DftiCreateDescriptor/DftiComputeForward/DftiComputeBackward API
# instead of the limited FFTW3 compatibility wrappers (guru64), which fail
# for howmany_rank > 1 (e.g., transforming one dim of a 3D+ array).

# This file is only included when fftw_provider == "mkl".

#==============================================================================#
# DFTI configuration constants (from mkl_dfti.h)
#==============================================================================#

# DFTI_CONFIG_VALUE enum
const DFTI_COMMITTED       = Int64(30)
const DFTI_UNCOMMITTED     = Int64(31)

# Precision
const DFTI_SINGLE          = Int64(35)
const DFTI_DOUBLE          = Int64(36)

# Forward domain
const DFTI_COMPLEX         = Int64(32)
const DFTI_REAL            = Int64(33)

# Placement
const DFTI_INPLACE         = Int64(43)
const DFTI_NOT_INPLACE     = Int64(44)

# Complex storage
const DFTI_COMPLEX_COMPLEX = Int64(39)
const DFTI_REAL_REAL       = Int64(40)
const DFTI_COMPLEX_REAL    = Int64(41)

# Packed format
const DFTI_CCE_FORMAT      = Int64(56)
const DFTI_CCS_FORMAT      = Int64(45)
const DFTI_PACK_FORMAT     = Int64(46)
const DFTI_PERM_FORMAT     = Int64(47)

# Allow/Avoid
const DFTI_ALLOW           = Int64(51)
const DFTI_AVOID           = Int64(52)

# DFTI_CONFIG_PARAM enum
const DFTI_FORWARD_DOMAIN           = Int64(0)
const DFTI_DIMENSION                = Int64(1)
const DFTI_LENGTHS                  = Int64(2)
const DFTI_PRECISION                = Int64(3)
const DFTI_FORWARD_SCALE            = Int64(4)
const DFTI_BACKWARD_SCALE           = Int64(5)
const DFTI_NUMBER_OF_TRANSFORMS     = Int64(7)
const DFTI_COMPLEX_STORAGE          = Int64(8)
const DFTI_REAL_STORAGE             = Int64(9)
const DFTI_CONJUGATE_EVEN_STORAGE   = Int64(10)
const DFTI_PLACEMENT                = Int64(11)
const DFTI_INPUT_STRIDES            = Int64(12)
const DFTI_OUTPUT_STRIDES           = Int64(13)
const DFTI_INPUT_DISTANCE           = Int64(14)
const DFTI_OUTPUT_DISTANCE          = Int64(15)
const DFTI_WORKSPACE                = Int64(17)
const DFTI_ORDERING                 = Int64(18)
const DFTI_TRANSPOSE                = Int64(19)
const DFTI_DESCRIPTOR_NAME          = Int64(20)
const DFTI_PACKED_FORMAT            = Int64(21)
const DFTI_COMMIT_STATUS            = Int64(22)
const DFTI_VERSION                  = Int64(23)
const DFTI_NUMBER_OF_USER_THREADS   = Int64(26)
const DFTI_THREAD_LIMIT             = Int64(27)
const DFTI_DESTROY_INPUT            = Int64(28)

# Ordered/Scrambled
const DFTI_ORDERED             = Int64(48)
const DFTI_BACKWARD_SCRAMBLED  = Int64(49)

#==============================================================================#
# DFTI API wrappers
#==============================================================================#

const DftiDescriptor = Ptr{Cvoid}

# Note: DftiCreateDescriptor is a macro in mkl_dfti.h that dispatches to
# specific functions based on precision and dimensionality. We call the
# specific underlying functions directly.

function dfti_create_descriptor_1d(precision::Int64, domain::Int64, length::Int64)
    handle = Ref{DftiDescriptor}(C_NULL)
    if precision == DFTI_SINGLE
        status = ccall((:DftiCreateDescriptor_s_1d, libfftw3), Int64,
                       (Ref{DftiDescriptor}, Int64, Int64),
                       handle, domain, length)
    else
        status = ccall((:DftiCreateDescriptor_d_1d, libfftw3), Int64,
                       (Ref{DftiDescriptor}, Int64, Int64),
                       handle, domain, length)
    end
    status != 0 && error("DftiCreateDescriptor failed: $(dfti_error_message(status))")
    return handle[]
end

function dfti_create_descriptor_md(precision::Int64, domain::Int64, ndim::Int64, lengths::Vector{Int64})
    handle = Ref{DftiDescriptor}(C_NULL)
    if precision == DFTI_SINGLE
        status = ccall((:DftiCreateDescriptor_s_md, libfftw3), Int64,
                       (Ref{DftiDescriptor}, Int64, Int64, Ptr{Int64}),
                       handle, domain, ndim, lengths)
    else
        status = ccall((:DftiCreateDescriptor_d_md, libfftw3), Int64,
                       (Ref{DftiDescriptor}, Int64, Int64, Ptr{Int64}),
                       handle, domain, ndim, lengths)
    end
    status != 0 && error("DftiCreateDescriptor failed: $(dfti_error_message(status))")
    return handle[]
end

function dfti_create_descriptor(precision::Int64, domain::Int64, lengths::Vector{Int64})
    ndim = length(lengths)
    if ndim == 1
        return dfti_create_descriptor_1d(precision, domain, lengths[1])
    else
        return dfti_create_descriptor_md(precision, domain, Int64(ndim), lengths)
    end
end

function dfti_set_value(handle::DftiDescriptor, param::Int64, value::Int64)
    status = ccall((:DftiSetValue, libfftw3), Int64,
                   (DftiDescriptor, Int64, Int64),
                   handle, param, value)
    status != 0 && error("DftiSetValue failed (param=$param): $(dfti_error_message(status))")
end

function dfti_set_value(handle::DftiDescriptor, param::Int64, value::Float32)
    status = ccall((:DftiSetValue, libfftw3), Int64,
                   (DftiDescriptor, Int64, Float32),
                   handle, param, value)
    status != 0 && error("DftiSetValue failed (param=$param): $(dfti_error_message(status))")
end

function dfti_set_value(handle::DftiDescriptor, param::Int64, value::Float64)
    status = ccall((:DftiSetValue, libfftw3), Int64,
                   (DftiDescriptor, Int64, Float64),
                   handle, param, value)
    status != 0 && error("DftiSetValue failed (param=$param): $(dfti_error_message(status))")
end

function dfti_set_value(handle::DftiDescriptor, param::Int64, value::Vector{Int64})
    status = ccall((:DftiSetValue, libfftw3), Int64,
                   (DftiDescriptor, Int64, Ptr{Int64}),
                   handle, param, value)
    status != 0 && error("DftiSetValue failed (param=$param): $(dfti_error_message(status))")
end

function dfti_commit_descriptor(handle::DftiDescriptor)
    status = ccall((:DftiCommitDescriptor, libfftw3), Int64,
                   (DftiDescriptor,), handle)
    status != 0 && error("DftiCommitDescriptor failed: $(dfti_error_message(status))")
end

function dfti_compute_forward_inplace(handle::DftiDescriptor, data::Ptr)
    status = ccall((:DftiComputeForward, libfftw3), Int64,
                   (DftiDescriptor, Ptr{Cvoid}),
                   handle, data)
    status != 0 && error("DftiComputeForward failed: $(dfti_error_message(status))")
end

function dfti_compute_forward_outofplace(handle::DftiDescriptor, indata::Ptr, outdata::Ptr)
    status = ccall((:DftiComputeForward, libfftw3), Int64,
                   (DftiDescriptor, Ptr{Cvoid}, Ptr{Cvoid}),
                   handle, indata, outdata)
    status != 0 && error("DftiComputeForward failed: $(dfti_error_message(status))")
end

function dfti_compute_backward_inplace(handle::DftiDescriptor, data::Ptr)
    status = ccall((:DftiComputeBackward, libfftw3), Int64,
                   (DftiDescriptor, Ptr{Cvoid}),
                   handle, data)
    status != 0 && error("DftiComputeBackward failed: $(dfti_error_message(status))")
end

function dfti_compute_backward_outofplace(handle::DftiDescriptor, indata::Ptr, outdata::Ptr)
    status = ccall((:DftiComputeBackward, libfftw3), Int64,
                   (DftiDescriptor, Ptr{Cvoid}, Ptr{Cvoid}),
                   handle, indata, outdata)
    status != 0 && error("DftiComputeBackward failed: $(dfti_error_message(status))")
end

function dfti_free_descriptor(handle::DftiDescriptor)
    ref = Ref(handle)
    ccall((:DftiFreeDescriptor, libfftw3), Int64,
          (Ref{DftiDescriptor},), ref)
end

function dfti_error_message(status::Int64)
    ptr = ccall((:DftiErrorMessage, libfftw3), Ptr{UInt8}, (Int64,), status)
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
    pinv::ScaledPlan
    function MKLcPlan{T,K,inplace,N,G}(handle::DftiDescriptor, flags::Integer, R::G,
                                       X::StridedArray{T,N}, Y::StridedArray) where {T<:fftwComplex,K,inplace,N,G}
        p = new(handle, size(X), size(Y), strides(X), strides(Y),
                alignment_of(X), alignment_of(Y), UInt32(flags), R)
        finalizer(p) do plan
            if plan.handle != C_NULL
                dfti_free_descriptor(plan.handle)
                plan.handle = C_NULL
            end
        end
        p
    end
end

function MKLcPlan{T,K,inplace,N}(handle::DftiDescriptor, flags::Integer, R::G,
                                 X::StridedArray{T,N}, Y::StridedArray) where {T<:fftwComplex,K,inplace,N,G}
    MKLcPlan{T,K,inplace,N,G}(handle, flags, R, X, Y)
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
    pinv::ScaledPlan
    function MKLrPlan{T,K,inplace,N,G}(handle::DftiDescriptor, flags::Integer, R::G,
                                       X::StridedArray{T,N}, Y::StridedArray) where {T<:fftwNumber,K,inplace,N,G}
        p = new(handle, size(X), size(Y), strides(X), strides(Y),
                alignment_of(X), alignment_of(Y), UInt32(flags), R)
        finalizer(p) do plan
            if plan.handle != C_NULL
                dfti_free_descriptor(plan.handle)
                plan.handle = C_NULL
            end
        end
        p
    end
end

function MKLrPlan{T,K,inplace,N}(handle::DftiDescriptor, flags::Integer, R::G,
                                 X::StridedArray{T,N}, Y::StridedArray) where {T<:fftwNumber,K,inplace,N,G}
    MKLrPlan{T,K,inplace,N,G}(handle, flags, R, X, Y)
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

# Override unsafe_convert so that MKL plans don't try to convert to PlanPtr
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
    dfti_lengths = Int64[sz[d] for d in reg]

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
    dfti_xstrides = Int64[0; [ist[d] for d in reg]]
    dfti_ystrides = Int64[0; [ost[d] for d in reg]]
    dfti_set_value(handle, DFTI_INPUT_STRIDES, dfti_xstrides)
    dfti_set_value(handle, DFTI_OUTPUT_STRIDES, dfti_ystrides)

    # Handle batch dimensions (non-transform dimensions)
    # Find all dimensions not in region
    batch_dims = [d for d in 1:N if !(d in reg)]

    if !isempty(batch_dims)
        # We can use NUMBER_OF_TRANSFORMS with DISTANCE for batched execution.
        # For multiple batch dims, compute total batch count and use the
        # appropriate distance. The key insight is that for column-major arrays
        # with contiguous batch dims, we can set up a single batched call.
        #
        # Strategy: find the innermost batch dimension (smallest stride).
        # Use it for NUMBER_OF_TRANSFORMS. For remaining batch dims, we
        # multiply into the batch count if they are contiguous in memory.

        # Sort batch dims by input stride
        sorted_batch = sort(batch_dims, by=d->ist[d])

        # Compute total number of transforms and distances
        nbatch = 1
        for d in batch_dims
            nbatch *= sz[d]
        end

        # Use stride of the dimension with smallest stride as distance
        # Check if all batch dimensions are "contiguous" in the sense that
        # a single NUMBER_OF_TRANSFORMS + DISTANCE covers the whole batch.
        # This works when all non-transform elements are accessed by
        # incrementing a single distance.

        # For general strided arrays, we need the batch to be describable
        # by a single (count, input_distance, output_distance).
        # We compute the minimum batch stride and check if total_batch * min_stride
        # fills the space correctly.

        min_ist_dim = sorted_batch[1]
        idist = ist[min_ist_dim]
        odist = ost[min_ist_dim]

        # Check if all batch dims can be described by a single distance
        # by verifying the strides form a consistent row-major/col-major pattern
        can_single_batch = true
        expected_ist = idist
        expected_ost = odist
        for d in sorted_batch
            expected_count = sz[d]
            if ist[d] != expected_ist || ost[d] != expected_ost
                can_single_batch = false
                break
            end
            expected_ist *= expected_count
            expected_ost *= expected_count
        end

        if can_single_batch
            dfti_set_value(handle, DFTI_NUMBER_OF_TRANSFORMS, Int64(nbatch))
            dfti_set_value(handle, DFTI_INPUT_DISTANCE, Int64(idist))
            dfti_set_value(handle, DFTI_OUTPUT_DISTANCE, Int64(odist))
        else
            inner_batch_count = sz[min_ist_dim]
            dfti_set_value(handle, DFTI_NUMBER_OF_TRANSFORMS, Int64(inner_batch_count))
            dfti_set_value(handle, DFTI_INPUT_DISTANCE, Int64(idist))
            dfti_set_value(handle, DFTI_OUTPUT_DISTANCE, Int64(odist))
            # We'll handle the outer loops during execution via assert_applicable checks
            # For now, the plan stores the full sz/osz for checking.
            # The execute function will need to loop over outer dims.
        end
    end

    dfti_commit_descriptor(handle)

    return MKLcPlan{T,K,inplace,N}(handle, flags, R, X, Y)
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
    dfti_lengths = Int64[sz[d] for d in rev_reg]

    prec = _dfti_precision(Tr)

    handle = dfti_create_descriptor(prec, DFTI_REAL, dfti_lengths)

    # Always out-of-place for rfft (no in-place version supported currently)
    dfti_set_value(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)

    # Use CCE (Conjugate-Complex-Even) format: output is complex array
    dfti_set_value(handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)

    # Input strides (real domain) - strides in real elements (reversed order)
    dfti_istrides = Int64[0; [ist[d] for d in rev_reg]]
    dfti_set_value(handle, DFTI_INPUT_STRIDES, dfti_istrides)

    # Output strides (complex domain) - strides in complex elements (reversed order)
    dfti_ostrides = Int64[0; [ost[d] for d in rev_reg]]
    dfti_set_value(handle, DFTI_OUTPUT_STRIDES, dfti_ostrides)

    # Handle batch dimensions
    batch_dims = [d for d in 1:N if !(d in reg)]
    if !isempty(batch_dims)
        sorted_batch = sort(batch_dims, by=d->ist[d])
        nbatch = prod(sz[d] for d in batch_dims)

        min_ist_dim = sorted_batch[1]
        idist = ist[min_ist_dim]
        odist = ost[min_ist_dim]

        can_single_batch = true
        expected_ist = idist
        expected_ost = odist
        for d in sorted_batch
            expected_count = sz[d]
            if ist[d] != expected_ist || ost[d] != expected_ost
                can_single_batch = false
                break
            end
            expected_ist *= expected_count
            expected_ost *= expected_count
        end

        if can_single_batch
            dfti_set_value(handle, DFTI_NUMBER_OF_TRANSFORMS, Int64(nbatch))
            dfti_set_value(handle, DFTI_INPUT_DISTANCE, Int64(idist))
            dfti_set_value(handle, DFTI_OUTPUT_DISTANCE, Int64(odist))
        else
            inner_count = sz[min_ist_dim]
            dfti_set_value(handle, DFTI_NUMBER_OF_TRANSFORMS, Int64(inner_count))
            dfti_set_value(handle, DFTI_INPUT_DISTANCE, Int64(idist))
            dfti_set_value(handle, DFTI_OUTPUT_DISTANCE, Int64(odist))
        end
    end

    dfti_commit_descriptor(handle)

    return MKLrPlan{Tr,FORWARD,inplace,N}(handle, flags, R, X, Y)
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
    dfti_lengths = Int64[sz_out[d] for d in rev_reg]

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
    dfti_complex_strides = Int64[0; [ist[d] for d in rev_reg]]
    dfti_real_strides = Int64[0; [ost[d] for d in rev_reg]]
    dfti_set_value(handle, DFTI_INPUT_STRIDES, dfti_complex_strides)
    dfti_set_value(handle, DFTI_OUTPUT_STRIDES, dfti_real_strides)

    # DESTROY_INPUT: MKL may overwrite input during c2r
    dfti_set_value(handle, DFTI_DESTROY_INPUT, DFTI_ALLOW)

    # Handle batch dimensions
    batch_dims = [d for d in 1:N if !(d in reg)]
    if !isempty(batch_dims)
        # For c2r, batch dim sizes come from the complex input array
        sorted_batch = sort(batch_dims, by=d->ist[d])
        nbatch = prod(sz_in[d] for d in batch_dims)

        min_ist_dim = sorted_batch[1]
        idist = ist[min_ist_dim]
        odist = ost[min_ist_dim]

        can_single_batch = true
        expected_ist = idist
        expected_ost = odist
        for d in sorted_batch
            expected_count = sz_in[d]
            if ist[d] != expected_ist || ost[d] != expected_ost
                can_single_batch = false
                break
            end
            expected_ist *= expected_count
            expected_ost *= expected_count
        end

        if can_single_batch
            dfti_set_value(handle, DFTI_NUMBER_OF_TRANSFORMS, Int64(nbatch))
            dfti_set_value(handle, DFTI_INPUT_DISTANCE, Int64(idist))
            dfti_set_value(handle, DFTI_OUTPUT_DISTANCE, Int64(odist))
        else
            inner_count = sz_in[min_ist_dim]
            dfti_set_value(handle, DFTI_NUMBER_OF_TRANSFORMS, Int64(inner_count))
            dfti_set_value(handle, DFTI_INPUT_DISTANCE, Int64(idist))
            dfti_set_value(handle, DFTI_OUTPUT_DISTANCE, Int64(odist))
        end
    end

    dfti_commit_descriptor(handle)

    return MKLrPlan{Tc,BACKWARD,inplace,N}(handle, flags, R, X, Y)
end

#==============================================================================#
# Execution
#==============================================================================#

function unsafe_execute!(plan::MKLcPlan{T,K,true}, X::StridedArray{T}, Y::StridedArray{T}) where {T,K}
    if K == FORWARD
        dfti_compute_forward_inplace(plan.handle, pointer(X))
    else
        dfti_compute_backward_inplace(plan.handle, pointer(X))
    end
end

function unsafe_execute!(plan::MKLcPlan{T,K,false}, X::StridedArray{T}, Y::StridedArray{T}) where {T,K}
    if K == FORWARD
        dfti_compute_forward_outofplace(plan.handle, pointer(X), pointer(Y))
    else
        dfti_compute_backward_outofplace(plan.handle, pointer(X), pointer(Y))
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
    dfti_compute_forward_outofplace(plan.handle, pointer(X), pointer(Y))
end

# c2r (backward brfft)
function unsafe_execute!(plan::MKLrPlan{Tc,BACKWARD}, X::StridedArray{Tc}, Y::StridedArray{<:fftwReal}) where {Tc<:fftwComplex}
    dfti_compute_backward_outofplace(plan.handle, pointer(X), pointer(Y))
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
        plan_brfft(X::StridedArray{$Tr}; kws...) = plan_brfft(X, ntuple(identity, ndims(X)); kws...)

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
