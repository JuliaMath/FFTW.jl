if VERSION < v"0.7.0-DEV.602"
    import Base.FFTW
else
    include("package.jl") # defines module FFTW
end
