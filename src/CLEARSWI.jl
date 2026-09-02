module CLEARSWI

using MriResearchTools
using Statistics

# Evaluated while this package is precompiled, so the version is part of the
# image and does not depend on path metadata being readable at run time. See
# package_version, which prefers it over pkgversion.
const PKG_VERSION = pkgversion(@__MODULE__)

include("tissue.jl")
include("utility.jl")
include("functions.jl")
include("magnitude_processing.jl")
include("phase_processing.jl")
include("simulate_single_echo.jl")

# The reference for the method this package implements, registered next to the
# code that implements it. The writer and the registry live in ROMEO, which has
# no dependencies, so nothing here points back up the dependency graph.
function __init__()
    register_citation!(:clearswi,
        """Eckstein, K., Bachrata, B., Hangel, G., Widhalm, G., Enzinger, C., Barth, M., Trattnig, S., Robinson, S.D., 2021.
           Improved susceptibility weighted imaging at ultra-high field using bipolar multi-echo acquisition and optimized image processing: CLEAR-SWI.
           NeuroImage 237, 118175.
           https://doi.org/10.1016/j.neuroimage.2021.118175""";
        label = "CLEAR-SWI")
end

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
