import Pkg
Pkg.activate(@__DIR__)
try
    using CLEARSWI, QuantitativeSusceptibilityMappingTGV, ArgParse
catch
    Pkg.add(["CLEARSWI", "QuantitativeSusceptibilityMappingTGV", "ArgParse"])
    using CLEARSWI, QuantitativeSusceptibilityMappingTGV, ArgParse
end

@time msg = clearswi_main(ARGS)
println(msg)
