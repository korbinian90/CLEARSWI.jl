using CLEARSWI
using Test
using TestItemRunner

@testset "CLEARSWI.jl" begin 
    @testset "Utils Tests" begin include("utility_test.jl") end
    @testset "Functions Test" begin include("functions_test.jl") end
    @testset "With FFTW Test" begin include("fftw_test.jl") end
end

@testset "ClearswiApp.jl" begin
    using ArgParse
    include("ClearswiApp.jl")
end

@testset "QSM" begin
    include("qsm.jl")
end

## print version to verify
println()
clearswi_main(["--version"])
