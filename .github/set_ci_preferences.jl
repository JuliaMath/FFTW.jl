#!/usr/bin/env julia 

open(joinpath(dirname(@__DIR__), "test", "Project.toml"), "a") do io
    println(io, """
    [preferences.FFTW]
    provider = "$(ARGS[1])"
    """)
end