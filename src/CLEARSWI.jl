
module CLEARSWI

using MriResearchTools
using Statistics

include("tissue.jl")
include("utility.jl")
include("functions.jl")

export calculateSWI,
        createMIP,
        saveconfiguration,
        Data,
        Options,
        savenii,
        readmag,
        readphase

"""
    CLEARSWI.dir(path...)
Construct a path relative to SWI root.
# Example
```julia
julia> CLEARSWI.dir("test","testData","small","Mag.nii")
"/home/korbinian90/.julia/dev/CLEARSWI/test/testData/small/Mag.nii"
```
"""
dir(path...) = joinpath(dirname(@__DIR__),path...)

end
