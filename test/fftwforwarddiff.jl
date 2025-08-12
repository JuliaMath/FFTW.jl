using FFTW, ForwardDiff, Test
using ForwardDiff: Dual, value, partials

@testset "ForwardDiff extension" begin
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
        @test ForwardDiff.derivative(f, 0.1) ≡ 1.0

        @test mul!(similar(x1), FFTW.plan_r2r(x1, FFTW.R2HC), x1) == FFTW.r2r(x1, FFTW.R2HC)

        x = [Dual(1.0,2,3), Dual(4,5,6)]
        a = FFTW.r2r(x, FFTW.REDFT00)
        b = FFTW.r2r!(x, FFTW.REDFT00)
        @test a == b == x
    end

    @testset "dct" begin
        x = [Dual(1.0,2,3), Dual(4,5,6)]
        a = dct(x)
        b = dct!(x)
        @test a == b == x

        c = x -> dct([x; 0; 0])[1]
        @test ForwardDiff.derivative(c,0.1) ≈ 1/sqrt(3)
    end
end