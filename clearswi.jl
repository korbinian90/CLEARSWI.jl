import Pkg
Pkg.activate(@__DIR__)
using CLEARSWI, ArgParse

@time msg = clearswi_main(ARGS)
println(msg)
