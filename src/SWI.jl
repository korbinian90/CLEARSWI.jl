
module SWI

using MRI, Statistics

include("tissue.jl")
include("utility.jl")
include("functions.jl")

export calculateSWI, createMIP, saveconfiguration, Data, Options

"""
    SWI.dir(path...)
Construct a path relative to SWI root.
# Example
```julia
julia> SWI.dir("test","testData","small","Mag.nii")
"/home/korbinian90/.julia/dev/SWI/test/testData/small/Mag.nii"
```
"""
dir(path...) = joinpath(dirname(@__DIR__),path...)

end
