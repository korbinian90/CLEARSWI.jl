[![Build Status](https://travis-ci.com/korbinian90/SWI.jl.svg?branch=master)](https://travis-ci.com/korbinian90/SWI.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/korbinian90/SWI.jl?svg=true)](https://ci.appveyor.com/project/korbinian90/SWI-jl)
[![Codecov](https://codecov.io/gh/korbinian90/SWI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/SWI.jl)
[![Coveralls](https://coveralls.io/repos/github/korbinian90/SWI.jl/badge.svg?branch=master)](https://coveralls.io/github/korbinian90/SWI.jl?branch=master)

# Susceptibility Weighted Imaging (SWI)
SWI provides magnetic resonance images with improved vein and iron contrast by weighting the magnitude image with a preprocessed phase image. This package has the additional capability of multi-echo SWI, intensity correction and improved phase processing. It can also work with classical single-echo SWI.

## Getting Started

### Prerequisites
A Julia installation v1.x

Magnitude and Phase images in NIfTI fileformat (4D images with echoes in the 4th dimension)

### Installing
Open the REPL in Julia and type

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/korbinian90/MriResearchTools.jl"))
Pkg.add(PackageSpec(url="https://github.com/korbinian90/SWI.jl"))
```

### Usage
This is a simple multi-echo case without changing default behavior
```julia
using SWI

TEs = [4,8,12] # change this to the Echo Time of your SWI sequence. For multi-echoes, set a list of TE values, else set a list with a single TE value.
nifti_folder = SWI.dir("test","testData","small")
magfile = joinpath(nifti_folder, "Mag.nii") # Path to the magnitude image in nifti format, must be .nii or .hdr
phasefile = joinpath(nifti_folder, "Phase.nii") # Path to the phase image

mag = readmag(magfile)
phase = readphase(phasefile)
data = Data(mag, phase, mag.header, TEs)

swi = calculateSWI(data)
mip = createMIP(swi)

savenii(swi, "<outputpath>/swi.nii"; header=mag.header) # change <outputpath> with the path where you want to save the reconstructed SWI images
savenii(mip, "<outputpath>/mip.nii"; header=mag.header)
```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/SWI.jl/blob/development/LICENSE) for details
