function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    for n = 1:3, T in (Float32, Float64, ComplexF32, ComplexF64), D in (UnitRange{Int}, Vector{Int}, Int)
        precompile(Tuple{typeof(fft),Array{T,n},D})
        precompile(Tuple{typeof(rfft),Array{T,n},D})
        precompile(Tuple{typeof(ifft),Array{T,n},D})
        precompile(Tuple{typeof(irfft),Array{T,n},Int,D})
        precompile(Tuple{typeof(rfft_output_size),Tuple{Int, Int, Int},D})
        precompile(Tuple{typeof(rfft_output_size),Tuple{Int, Int},D})
        precompile(Tuple{typeof(rfft_output_size),Tuple{Int},D})
    end
end
