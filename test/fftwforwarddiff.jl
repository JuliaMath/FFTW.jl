@testset "r2r" begin
    x1 = Dual.(1:4.0, 2:5, 3:6)
    t = FFTW.r2r(x1, FFTW.R2HC)

    @test value.(t) == FFTW.r2r(value.(x1), FFTW.R2HC)
    @test partials.(t, 1) == FFTW.r2r(partials.(x1, 1), FFTW.R2HC)
    @test partials.(t, 2) == FFTW.r2r(partials.(x1, 2), FFTW.R2HC)

    t = FFTW.r2r(x1 + 2im*x1, FFTW.R2HC)
    @test value.(t) == FFTW.r2r(value.(x1 + 2im*x1), FFTW.R2HC)
    @test partials.(t, 1) == FFTW.r2r(partials.(x1 + 2im*x1, 1), FFTW.R2HC)
    @test partials.(t, 2) == FFTW.r2r(partials.(x1 + 2im*x1, 2), FFTW.R2HC)

    f = ω -> FFTW.r2r([ω; zeros(9)], FFTW.R2HC)[1]
    @test derivative(f, 0.1) ≡ 1.0

    @test mul!(similar(x1), FFTW.plan_r2r(x1, FFTW.R2HC), x1) == FFTW.r2r(x1, FFTW.R2HC)
end