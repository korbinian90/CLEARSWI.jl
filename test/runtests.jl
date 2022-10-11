using CLEARSWI
using Test
using TestItemRunner

@time begin
    @time @testset "Utils Tests" begin include("utility_test.jl") end
    @time @testset "Functions Test" begin include("functions_test.jl") end
    @time @testset "With FFTW Test" begin include("fftw_test.jl") end
    include("qsm.jl")
end
