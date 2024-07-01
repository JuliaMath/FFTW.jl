let a = rand(Float64,(8,4,4)), b = PaddedRFFTArray(a), c = copy(b)

@testset "PaddedRFFTArray creation" begin
  @test a == real(b)
  @test c == b
  @test c.r == b.r
  @test typeof(similar(b)) === typeof(b)
  @test size(similar(b,Float32)) === size(b)
  @test size(similar(b,Float32).r) === size(b.r)
  @test size(similar(b,(4,4,4)).r) === (4,4,4)
  @test size(similar(b,Float32,(4,4,4)).r) === (4,4,4) 
end

@testset "rfft! and irfft!" begin
  @test rfft(a) ≈ rfft!(b) 
  @test a ≈ irfft!(b)
  @test rfft(a,1:2) ≈ rfft!(b,1:2) 
  @test a ≈ irfft!(b,1:2)
  @test rfft(a,(1,3)) ≈ rfft!(b,(1,3)) 
  @test a ≈ irfft!(b,(1,3))
  
  p = plan_rfft!(c)
  @test p*c ≈ rfft!(b)
  @test p\c ≈ irfft!(b)

  a = rand(Float64,(9,4,4))
  b = PaddedRFFTArray(a)
  @test a == real(b)
  @test rfft(a) ≈ rfft!(b) 
  @test a ≈ irfft!(b)
  @test rfft(a,1:2) ≈ rfft!(b,1:2) 
  @test a ≈ irfft!(b,1:2)
  @test rfft(a,(1,3)) ≈ rfft!(b,(1,3)) 
  @test a ≈ irfft!(b,(1,3))
end

@testset "Read binary file to PaddedRFFTArray" begin
  for s in ((8,4,4),(9,4,4),(8,),(9,))
    aa = rand(Float64,s)
    f = IOBuffer()
    write(f,aa)
    @test aa == real(PaddedRFFTArray(seekstart(f),s))
    aa = rand(Float32,s)
    f = IOBuffer()
    write(f,aa)
    @test aa == real(PaddedRFFTArray{Float32}(seekstart(f),s))
  end
end

@testset "brfft!" begin
  a = rand(Float64,(4,4))
  b = PaddedRFFTArray(a)
  rfft!(b)
  @test (brfft!(b) ./ 16) ≈ a
end

@testset "FFTW MEASURE flag" begin
  c = similar(b)
  p = plan_rfft!(c,flags=FFTW.MEASURE)
  p.pinv = plan_irfft!(c,flags=FFTW.MEASURE)
  c .= b 
  @test c == b
  @test p*c ≈ rfft!(b)
  @test p\c ≈ irfft!(b)
end
end #let block
