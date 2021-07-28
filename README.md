[![Build Status](https://github.com/korbinian90/CLEARSWI.jl/workflows/CI/badge.svg)](https://github.com/korbinian90/CLEARSWI.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/korbinian90/CLEARSWI.jl?svg=true)](https://ci.appveyor.com/project/korbinian90/CLEARSWI-jl)
[![Codecov](https://codecov.io/gh/korbinian90/CLEARSWI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/CLEARSWI.jl)

# Susceptibility Weighted Imaging (CLEAR-SWI)
Published at the [ISMRM as CLEAR-SWI](https://index.mirasmart.com/ISMRM2020/PDFfiles/3201.html). It provides magnetic resonance images with improved vein and iron contrast by weighting a combined magnitude image with a preprocessed phase image. This package has the additional capability of multi-echo SWI, intensity correction, contrast enhancement and improved phase processing. The reason for the development of this package was to solve artefacts at ultra-high field strength (7T), however, it also drastically improves the SWI quality at lower field strength.

## Getting Started

### Prerequisites
A Julia installation ≥ 1.3 ([Official Julia Webpage](https://julialang.org/downloads/))

Single-echo or multi-echo Magnitude and Phase images in NIfTI fileformat (4D images with echoes in the 4th dimension)

### Installing
Run the following commands in Julia (either interactively in the REPL or as a script)

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/korbinian90/CLEARSWI.jl"))
```

### Usage
This is a simple multi-echo case without changing default behavior
```julia
using CLEARSWI

TEs = [4,8,12] # change this to the Echo Time of your sequence. For multi-echoes, set a list of TE values, else set a list with a single TE value.
nifti_folder = CLEARSWI.dir("test","testData","small") # replace with path to your folder e.g. nifti_folder="/data/clearswi"
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
options = Options(phase_hp_σ=[10,10,5], mag_softplus=false)
swi = calculateSWI(data, options);
```

All the possible options with the default values are
```julia
mag_combine=:SNR
mag_sens=nothing
mag_softplus=true
phase_unwrap=:laplacian
phase_hp_σ=[4,4,0]
phase_scaling_type=:tanh
phase_scaling_strength=4
writesteps=nothing
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

* The `phase_hp_σ` is used for high-pass filtering and is given in voxel for the [x,y,z]-dimension.  

* `phase_scaling_type` is the scaling function to create the phase mask and can be `:tanh` for sigmoidal filtering, or `:positive`, `:negative`, and `:triangular` for traditional SWI application. If the scaling has the wrong sign, the phase input can be simply inverted by a minus sign: `phase = -readphase(phasefile);`

* `phase_scaling_strength` adjusts the strength of the created phase mask. A higher phase_scaling_strength is a stronger phase appearance. With a traditional SWI `phase_scaling_type` it corresponds to the power or number of phase mask multiplications.

* Set `writesteps` to the path, where intermediate steps should be saved, e.g. `writesteps="/tmp/clearswi_steps"`. If set to `nothing`, intermediate steps won't be saved.

### Calculating T2* and B0 maps on multi-echo datasets
T2* and B0 maps can be calculated using the package [MriResearchTools](https://github.com/korbinian90/MriResearchTools.jl):

#### Installing:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/korbinian90/MriResearchTools.jl"))
```

#### Usage:
With the previously defined variables `phase`, `mag` and `TEs`

```julia
using MriResearchTools

unwrapped = romeo(phase; mag=mag, TEs=TEs) # type ?romeo in REPL for options
B0 = calculateB0_unwrapped(unwrapped, mag, TEs) # inverse variance weighted

t2s = NumART2star(mag, TEs)
r2s = r2s_from_t2s(t2s)
```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/CLEARSWI.jl/blob/master/LICENSE) for details
