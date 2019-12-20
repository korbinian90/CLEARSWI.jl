
module SWI

using MRI, Statistics

include("tissue.jl")
include("utility.jl")
include("functions.jl")

export calculateSWI, createMIP, saveconfiguration, Data, Options

end
