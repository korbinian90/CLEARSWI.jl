[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://korbinian90.github.io/CLEARSWI.jl/dev)
[![Build Status](https://github.com/korbinian90/CLEARSWI.jl/workflows/CI/badge.svg)](https://github.com/korbinian90/CLEARSWI.jl/actions)
[![Codecov](https://codecov.io/gh/korbinian90/CLEARSWI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/CLEARSWI.jl)

![test_clear_swi_github](https://user-images.githubusercontent.com/1307522/194285019-60e0e0a3-1bf5-4563-86bd-4201de2be08b.png)

Apply CLEARSWI in the command line without Julia programming experience. This repository contains the CLEARSWI algorithm and a command line interface.

# Susceptibility Weighted Imaging (CLEAR-SWI)

Published as [CLEAR-SWI](https://doi.org/10.1016/j.neuroimage.2021.118175). It provides magnetic resonance images with improved vein and iron contrast by weighting a combined magnitude image with a preprocessed phase image. This package has the additional capability of multi-echo SWI, intensity correction, contrast enhancement and improved phase processing. The reason for the development of this package was to solve artefacts at ultra-high field strength (7T), however, it also drastically improves the SWI quality at lower field strength.

## Option 1: Download standalone executables

https://github.com/korbinian90/CompileMRI.jl/releases  
This package contains binaries with dependencies.  
For usage help, run the command line program in `bin/clearswi` without arguments or see [below](https://github.com/korbinian90/CLEARSWI.jl#command-line-help).

## Option 2: Usage - command line via Julia

Install Julia 1.9 or newer (https://julialang.org)  
Copy the file `clearswi.jl` from this repository to a convenient location. An alias for `clearswi` as `julia <path-to-file>/clearswi.jl` might be useful.

```bash
    $ julia <path-to-file>/clearswi.jl -p phase.nii -m mag.nii -t [2.1,4.2,6.3] -o results
```

On the first run, the dependencies will be installed automatically.

For an extended explanation of the command line interface see [below](https://github.com/korbinian90/CLEARSWI).

## Option 3: Usage - Julia

### Prerequisites

A Julia installation ≥ 1.6 ([Official Julia Webpage](https://julialang.org/downloads/))

Single-echo or multi-echo Magnitude and Phase images in NIfTI fileformat (4D images with echoes in the 4th dimension)

### Installing

Run the following commands in Julia (either interactively in the REPL or as a script)

```julia
import Pkg; Pkg.add("CLEARSWI")
```

### Function Reference

https://korbinian90.github.io/CLEARSWI.jl/dev

### Usage

This is a simple multi-echo case without changing default behavior

```julia
using CLEARSWI

TEs = [4,8,12] # change this to the Echo Time of your sequence. For multi-echoes, set a list of TE values, else set a list with a single TE value.
nifti_folder = CLEARSWI.dir("test","data","small") # replace with path to your folder e.g. nifti_folder="/data/clearswi"
magfile = joinpath(nifti_folder, "Mag.nii") # Path to the magnitude image in nifti format, must be .nii or .hdr
phasefile = joinpath(nifti_folder, "Phase.nii") # Path to the phase image

mag = readmag(magfile);
phase = readphase(phasefile);
data = Data(mag, phase, mag.header, TEs);

swi = calculateSWI(data);
# mip = createIntensityProjection(swi, minimum); # minimum intensity projection, other Julia functions can be used instead of minimum
mip = createMIP(swi); # shorthand for createIntensityProjection(swi, minimum)

savenii(swi, "<outputpath>/swi.nii"; header=mag.header) # change <outputpath> with the path where you want to save the reconstructed SWI
savenii(mip, "<outputpath>/mip.nii"; header=mag.header)
```

### Available Options

To apply custom options use the following keyword syntax (example to apply 3D high-pass filtering for the phase with the given kernel size and deactivate softplus magnitude scaling):

```julia
options = Options(phase_hp_sigma=[10,10,5], mag_softplus=false)
swi = calculateSWI(data, options);
```

All the possible options with the default values are

```julia
mag_combine=:SNR
mag_sens=nothing
mag_softplus=true
phase_unwrap=:laplacian
phase_hp_sigma=[4,4,0]
phase_scaling_type=:tanh
phase_scaling_strength=4
writesteps=nothing
qsm=false
```

* `mag_combine` selects the echo combination for the magnitude. Options are 
  * `:SNR`
  * `:average`
  * `:last` to select the last echo
  * `(:CNR => (:gm, :wm))` to optimize the contrast between two selected tissues with the possible tissues classes to select in `src\tissue.jl`, currently only working for 7T
  * `(:echo => 3)` to select the 3rd echo 
  * `(:closest => 15.3)` to select the echo that is closest to 15.3 ms 
  * `(:SE => 15.3)` to simulate the contrast that would be achieved using a corresponding single-echo scan with 15.3 ms echo time.

* If `mag_sens` is set to an array, it is used instead of CLEAR-SWI sensitivity estimation. This can also be set to `mag_sens=[1]` to use the constant sensitivity of 1 and effectively avoid sensitivity correction.

* To deactivate scaling of the combined magnitude with the softplus function, use `mag_softplus=false`.

* `phase_unwrap` is either `:laplacian`, `:romeo`, or `:laplacianslice` (slicewise laplacian unwrapping)

* The `phase_hp_sigma` is used for high-pass filtering and is given in voxel for the [x,y,z]-dimension.  

* `phase_scaling_type` is the scaling function to create the phase mask and can be `:tanh` or `:negativetanh` for sigmoidal filtering, or `:positive`, `:negative`, and `:triangular` for traditional SWI application.

* `phase_scaling_strength` adjusts the strength of the created phase mask. A higher phase_scaling_strength is a stronger phase appearance. With a traditional SWI `phase_scaling_type` it corresponds to the power or number of phase mask multiplications.

* Set `writesteps` to the path, where intermediate steps should be saved, e.g. `writesteps="/tmp/clearswi_steps"`. If set to `nothing`, intermediate steps won't be saved.

* [Experimental] Set `qsm` to true to use QSM processing for the phase contrast. This requires the additional use of a QSM package. Supported are either [`QSM`](https://github.com/kamesy/QSM.jl) or [`QuantitativeSusceptibilityMappingTGV`](https://github.com/korbinian90/QuantitativeSusceptibilityMappingTGV.jl). Before setting this option, you need load one of these packages with `using`.

### Calculating T2* and B0 maps on multi-echo datasets

T2* and B0 maps can be calculated using the package [MriResearchTools](https://github.com/korbinian90/MriResearchTools.jl):

#### Installing

B0 and R2*/T2* calculation requires the package `MriResearchTools`

```julia
using Pkg
Pkg.add("MriResearchTools")
```

#### Usage

With the previously defined variables `phase`, `mag` and `TEs`

```julia
using MriResearchTools

unwrapped = romeo(phase; mag=mag, TEs=TEs) # type ?romeo in REPL for options
B0 = calculateB0_unwrapped(unwrapped, mag, TEs) # inverse variance weighted

t2s = NumART2star(mag, TEs)
r2s = r2s_from_t2s(t2s)
```

## Command Line Help

```plain
$ .\bin\clearswi
usage: <PROGRAM> [-m MAGNITUDE] [-p PHASE] [-o OUTPUT]
                 [-t ECHO-TIMES [ECHO-TIMES...]] [--qsm]
                 [--mag-combine MAG-COMBINE [MAG-COMBINE...]]
                 [--mag-sensitivity-correction MAG-SENSITIVITY-CORRECTION]
                 [--mag-softplus-scaling MAG-SOFTPLUS-SCALING]
                 [--unwrapping-algorithm UNWRAPPING-ALGORITHM]
                 [--filter-size FILTER-SIZE [FILTER-SIZE...]]
                 [--phase-scaling-type PHASE-SCALING-TYPE]
                 [--phase-scaling-strength PHASE-SCALING-STRENGTH]
                 [-e ECHOES [ECHOES...]] [-N] [--no-phase-rescale]
                 [--writesteps WRITESTEPS] [-v] [--version] [-h]

4.0.4

optional arguments:
  -m, --magnitude MAGNITUDE
                        The magnitude image (single or multi-echo)
  -p, --phase PHASE     The phase image (single or multi-echo)
  -o, --output OUTPUT   The output path or filename (default:
                        "clearswi.nii")
  -t, --echo-times ECHO-TIMES [ECHO-TIMES...]
                        The echo times are required for multi-echo
                        datasets specified in array or range syntax
                        (eg. "[1.5,3.0]" or "3.5:3.5:14").
  --qsm                 When activated uses RTS QSM for phase weighting
  --mag-combine MAG-COMBINE [MAG-COMBINE...]
                        SNR | average | echo <n> | SE <te>. Magnitude
                        combination algorithm. echo <n> selects a
                        specific echo; SE <te> simulates a single echo
                        scan of the given echo time. (default:
                        ["SNR"])
  --mag-sensitivity-correction MAG-SENSITIVITY-CORRECTION
                        <filename> | on | off. Use the CLEAR-SWI
                        sensitivity correction. Alternatively, a
                        sensitivity map can be read from a file
                        (default: "on")
  --mag-softplus-scaling MAG-SOFTPLUS-SCALING
                        on | off. Set softplus scaling of the
                        magnitude (default: "on")
  --unwrapping-algorithm UNWRAPPING-ALGORITHM
                        laplacian | romeo | laplacianslice (default:
                        "laplacian")
  --filter-size FILTER-SIZE [FILTER-SIZE...]
                        Size for the high-pass phase filter in voxels.
                        Can be given as <x> <y> <z> or in array syntax
                        (e.g. [2.2,3.1,0], which is effectively a 2D
                        filter). (default: ["[4,4,0]"])
  --phase-scaling-type PHASE-SCALING-TYPE
                        tanh | negativetanh | positive | negative |
                        triangular Select the type of phase scaling.
                        positive or negative with a strength of 3-6 is
                        used in standard SWI. (default: "tanh")
  --phase-scaling-strength PHASE-SCALING-STRENGTH
                        Sets the phase scaling strength. Corresponds
                        to power values for positive, negative and
                        triangular phase scaling type. (default: "4")
  -e, --echoes ECHOES [ECHOES...]
                        Load only the specified echoes from disk
                        (default: [":"])
  -N, --no-mmap         Deactivate memory mapping. Memory mapping
                        might cause problems on network storage
  --no-phase-rescale    Deactivate automatic rescaling of phase
                        images. By default the input phase is rescaled
                        to the range [-π;π].
  --writesteps WRITESTEPS
                        Set to the path of a folder, if intermediate
                        steps should be saved.
  -v, --verbose         verbose output messages
  --version             show version information and exit
  -h, --help            show this help message and exit
```

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/CLEARSWI.jl/blob/master/LICENSE) for details
