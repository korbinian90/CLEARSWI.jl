using SWI
using Test

@time begin
    @time @testset "Utils Tests" begin include("utility_test.jl") end
    @time @testset "Functions Test" begin include("functions_test.jl") end
end
