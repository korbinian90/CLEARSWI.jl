import Pkg
Pkg.activate(@__DIR__)
try
    using CLEARSWI, QuantitativeSusceptibilityMappingTGV
catch
    try
        Pkg.add("CLEARSWI")
    catch LoadError
        println("Skipping CLEARSWI installation, probably local directory used.")
    end
    Pkg.add("QuantitativeSusceptibilityMappingTGV")
    using CLEARSWI, QuantitativeSusceptibilityMappingTGV
end

@time msg = clearswi_main(ARGS)
println(msg)
