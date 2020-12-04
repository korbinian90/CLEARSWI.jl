module CLEARSWI

using MriResearchTools
using Statistics
using ImageFiltering
using LsqFit
using DSP

include("tissue.jl")
include("utility.jl")
include("functions.jl")
include("GEPCI.jl")
include("QuinnSWI.jl")

export calculateSWI,
        createMIP,
        saveconfiguration,
        Data,
        Options,
        savenii,
        readmag,
        readphase,
        GEPCI, quinn

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
