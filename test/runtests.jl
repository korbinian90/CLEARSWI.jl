using CLEARSWI
using TestItemRunner

@run_package_tests

@testitem "CLEARSWI.jl" begin 
    using Statistics
    @testset "Utils Tests" begin include("utility_test.jl") end
    @testset "Functions Test" begin include("functions_test.jl") end
    @testset "With FFTW Test" begin include("fftw_test.jl") end
end

@testitem "ClearswiApp.jl" begin
    include("ClearswiApp.jl")
end

@testitem "QSM" begin
    include("qsm.jl")
end

## print version to verify
println()
clearswi_main(["--version"])
