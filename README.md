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
Pkg.add(PackageSpec(url="https://github.com/korbinian90/CLEARSWI.jl", rev="other_me_methods"))
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

options = Options() # default CLEAR-SWI
#options = Options(combination_type=:GEPCI)
#options = Options(combination_type=:QuinnSWI)

swi = calculateSWI(data, options)
mip = createMIP(swi)

savenii(swi, "<outputpath>/swi.nii"; header=mag.header) # change <outputpath> with the path where you want to save the reconstructed SWI
savenii(mip, "<outputpath>/mip.nii"; header=mag.header)
```

### Available Options
The default options are
```julia
Options(;σ=[4,4,0], unwrapping=:laplacian, mode=:tanh, level=4, combination_type=:SNR, sensitivity=nothing, writedir=nothing, writesteps=false, magscale=identity)
```
* The `σ` is used for high-pass filtering and is given in voxel for the [x,y,z]-dimension.  
`unwrapping` is either `:laplacian`, `:romeo`, or `:laplacianslice` (slicewise laplacian unwrapping)

* `mode` is the scaling function to create the phase mask and can be `:tanh` for sigmoidal filtering, or `:positive`, `:negative`, and `:triangular` for traditional SWI application.

* `level` adjusts the strength of the created phase mask. A higher level is a stronger phase appearance. With a traditional SWI `mode` it corresponds to the power or number of phase mask multiplications.

* `combination_type` selects the echo combination for the magnitude. Options are `:SNR`; `:average`; `:last` to select the last echo; `(:CNR => (:gm, :wm))` to optimize the contrast between two selected tissues with the possible tissues classes to select in `src\tissue.jl`; `(:echo => 3)` to select the 3rd echo; `(:closest => 15.3)` to select the echo that is closest to 15.3 ms; `(:SE => 15.3)` to simulate the contrast that would be achieved using a corresponding single-echo scan with 15.3 ms echo time.

* Set `writesteps` to `true`, if you want the intermediate steps to be saved in the folder specified with `writedir="C:/tmp/clearswi_steps"`.

* If `sensitivity` is set to an array, it is used instead of CLEAR-SWI sensitivity estimation. This can also be set to `sensitvity=[1]` to use the constant sensitivity of 1 and effectively avoid sensitivity correction.

* `magscale` makes it possible to use an additional scaling function after magnitude echo combination. In the publication, the following softplus function was definedand was used as `magscale=spq2`:
```julia
function spq2(X)
    q = estimatequantile(X, 0.8)
    CLEARSWI.softplus.(X, q/2, 2)
end
```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/CLEARSWI.jl/blob/development/LICENSE) for details
