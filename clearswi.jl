import Pkg
Pkg.activate(@__DIR__)
try
    using CLEARSWI, QuantitativeSusceptibilityMappingTGV, ArgParse
catch
    try
        Pkg.add("CLEARSWI")
    catch LoadError
        println("Skipping CLEARSWI installation, probably local directory used.")
    end
    Pkg.add(["QuantitativeSusceptibilityMappingTGV", "ArgParse"])
    using CLEARSWI, QuantitativeSusceptibilityMappingTGV, ArgParse
end

@time msg = clearswi_main(ARGS)
println(msg)
