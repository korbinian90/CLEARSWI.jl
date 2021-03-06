module CLEARSWI

using MriResearchTools
using Statistics

include("tissue.jl")
include("utility.jl")
include("functions.jl")
include("magnitude_processing.jl")
include("phase_processing.jl")
include("simulate_single_echo.jl")

export calculateSWI,
        createMIP,
        createIntensityProjection,
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
