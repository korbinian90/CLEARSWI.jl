[![Build Status](https://travis-ci.com/korbinian90/CLEARSWI.jl.svg?branch=master)](https://travis-ci.com/korbinian90/CLEARSWI.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/korbinian90/CLEARSWI.jl?svg=true)](https://ci.appveyor.com/project/korbinian90/CLEARSWI-jl)
[![Codecov](https://codecov.io/gh/korbinian90/CLEARSWI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/CLEARSWI.jl)
[![Coveralls](https://coveralls.io/repos/github/korbinian90/CLEARSWI.jl/badge.svg?branch=master)](https://coveralls.io/github/korbinian90/CLEARSWI.jl?branch=master)

# Susceptibility Weighted Imaging (CLEAR SWI)
Published at the [ISMRM as CLEAR SWI](https://index.mirasmart.com/ISMRM2020/PDFfiles/3201.html). It provides magnetic resonance images with improved vein and iron contrast by weighting the magnitude image with a preprocessed phase image. This package has the additional capability of multi-echo SWI, intensity correction and improved phase processing. It can also work with classical single-echo SWI. It was developed to solve artefacts at ultra-high field strength (7T), however, it also drastically improves the SWI quality at lower field strength.

## This is the version used for r1 CLEAR-SWI publication
It includes GEPCI and Quinn et al. SWI used as multi-echo SWI comparison methods. GEPCI and Quinn et al. are devised for 3T application and are not perfectly suited for 7T application. At 7T Quinn et al. SWI is highly affected by homodyne filtering artefacts.  
Both Quinn et al. and CLEAR-SWI achieve higher SNR (and are faster) than GEPCI, because they use weighted combination using the known noise statistics instead of voxel-wise fitting. The fitting in GEPCI is additionally problematic as it assumes wrong minima in certain cases, although this might be avoided by a more careful selection of the fitting algorithm.

## Getting Started

### Prerequisites
A Julia installation ≥ 1.3

Magnitude and Phase images in NIfTI fileformat (4D images with echoes in the 4th dimension)

### Installing
Open the REPL in Julia and type

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/korbinian90/CLEARSWI.jl"))
```

### Usage
This is a simple multi-echo case without changing default behavior
```julia
using CLEARSWI

TEs = [4,8,12] # change this to the Echo Time of your sequence. For multi-echoes, set a list of TE values, else set a list with a single TE value.
nifti_folder = SWI.dir("test","testData","small")
magfile = joinpath(nifti_folder, "Mag.nii") # Path to the magnitude image in nifti format, must be .nii or .hdr
phasefile = joinpath(nifti_folder, "Phase.nii") # Path to the phase image

mag = readmag(magfile)
phase = readphase(phasefile)
data = Data(mag, phase, mag.header, TEs)

swi = calculateSWI(data)
mip = createMIP(swi)

savenii(swi, "<outputpath>/swi.nii"; header=mag.header) # change <outputpath> with the path where you want to save the reconstructed SWI
savenii(mip, "<outputpath>/mip.nii"; header=mag.header)
```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/CLEARSWI.jl/blob/development/LICENSE) for details
