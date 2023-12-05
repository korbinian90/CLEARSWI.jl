module CLEARSWI

using MriResearchTools
using Statistics

include("tissue.jl")
include("utility.jl")
include("functions.jl")
include("magnitude_processing.jl")
include("phase_processing.jl")
include("simulate_single_echo.jl")

clearswi_main(args...; kwargs...) = @warn("Type `using MriResearchTools ArgParse` to use this function \n `?clearswi_main` for argument help")
if !isdefined(Base, :get_extension) # fallback for julia < 1.9
    include("../ext/ClearswiApp/ClearswiApp.jl")
end

export calculateSWI,
    createMIP,
    createIntensityProjection,
    saveconfiguration,
    Data,
    Options,
    savenii,
    readmag,
    readphase,
    clearswi_main

"""
    CLEARSWI.dir(path...)
Construct a path relative to SWI root.
# Example
```julia
julia> CLEARSWI.dir("test","data","small","Mag.nii")
"/home/korbinian90/.julia/dev/CLEARSWI/test/data/small/Mag.nii"
```
"""
dir(path...) = joinpath(dirname(@__DIR__), path...)

end
