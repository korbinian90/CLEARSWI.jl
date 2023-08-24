import Pkg
Pkg.activate(@__DIR__)
try
    using CLEARSWI, ArgParse
catch
    Pkg.add(["CLEARSWI", "ArgParse"])
    using CLEARSWI, ArgParse
end

@time msg = clearswi_main(ARGS)
println(msg)
