import Pkg
Pkg.activate(@__DIR__)
try
    using CLEARSWI, QSM, ArgParse
catch
    Pkg.add(["CLEARSWI", "QSM", "ArgParse"])
    using CLEARSWI, QSM, ArgParse
end

@time msg = clearswi_main(ARGS)
println(msg)
