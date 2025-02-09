using ArgParse, QuantitativeSusceptibilityMappingTGV

original_path = abspath(".")
p = CLEARSWI.dir("test", "data", "small")
tmpdir = mktempdir()
cd(tmpdir)

phasefile_me = joinpath(p, "Phase.nii")
phasefile_me_nan = joinpath(p, "phase_with_nan.nii")
magfile_me = joinpath(p, "Mag.nii")
tmpdir = mktempdir()
phasefile_1eco = joinpath(tmpdir, "Phase.nii")
magfile_1eco = joinpath(tmpdir, "Mag.nii")
phasefile_1arreco = joinpath(tmpdir, "Phase.nii")
magfile_1arreco = joinpath(tmpdir, "Mag.nii")
magfile_me_nan_size = joinpath(tmpdir, "mag_nan_size.nii")
maskfile = joinpath(tmpdir, "Mask.nii")
savenii(readmag(magfile_me)[:,:,:,1] |> I -> I .> CLEARSWI.MriResearchTools.median(I), maskfile)
savenii(readphase(phasefile_me)[:,:,:,1], phasefile_1eco)
savenii(readmag(magfile_me)[:,:,:,1], magfile_1eco)
savenii(readphase(phasefile_me)[:,:,:,[1]], phasefile_1arreco)
savenii(readmag(magfile_me)[:,:,:,[1]], magfile_1arreco)
savenii(ones(size(readphase(phasefile_me_nan))), magfile_me_nan_size)

function test_clearswi(args)
    folder = tempname()
    args = [args..., "-o", folder]
    try
        msg = clearswi_main(args)
        @test msg == 0
        @test isfile(joinpath(folder, "clearswi.nii"))
    catch e
        println(args)
        println(sprint(showerror, e, catch_backtrace()))
        @test "test failed" == "with error" # signal a failed test
    end
end

configurations_se(pf, mf) = configurations_se(["-p", pf, "-m", mf])
configurations_se(pm) = [
    [pm...],
    # [pm..., "--qsm"],
    [pm..., "--mag-combine", "SNR"],
    [pm..., "--mag-combine", "average"],
    [pm..., "--mag-combine", "echo", "2"],
    [pm..., "--mag-combine", "SE", "2.4"],
    [pm..., "--mag-sensitivity-correction", "off"],
    [pm..., "--mag-softplus-scaling", "off"],
    [pm..., "--unwrapping-algorithm", "romeo"],
    #[pm..., "--unwrapping-algorithm", "laplacianslice"], # laplacianslice currently unstable
    [pm..., "--filter-size", "[2,2,3]"],
    [pm..., "--phase-scaling-type", "negativetanh"],
    [pm..., "--phase-scaling-type", "positive"],
    [pm..., "--phase-scaling-type", "negative"],
    [pm..., "--phase-scaling-type", "triangular"],
    [pm..., "--phase-scaling-strength", "1"],
    [pm..., "--phase-scaling-strength", "10"],
    [pm..., "-N"],
    [pm..., "--no-phase-rescale"],
    [pm..., "--writesteps", tmpdir],
    [pm..., "--mip-slices", "3"],
]
configurations_me(phasefile_me, magfile_me) = configurations_me(["-p", phasefile_me, "-m", magfile_me])
configurations_me(pm) = [
    [pm..., "-e", "1:2", "-t", "[2,4]"], # giving two echo times for two echoes used out of three
    [pm..., "-e", "[1,3]", "-t", "[2,4,6]"], # giving three echo times for two echoes used out of three
    [pm..., "-e", "[1", "3]", "-t", "[2,4,6]"],
    [pm..., "-t", "[2,4,6]"],
    [pm..., "-t", "2:2:6"],
    [pm..., "-t", "[2.1,4.2,6.3]"],
    # [pm..., "-t", "[2.1,4.2,6.3]", "--qsm"],
]

files = [(phasefile_1eco, magfile_1eco), (phasefile_1arreco, magfile_1arreco), (phasefile_1eco, magfile_1arreco), (phasefile_1arreco, magfile_1eco)]
for (pf, mf) in files, args in configurations_se(pf, mf)
    test_clearswi(args)
end
for args in configurations_me(phasefile_me, magfile_me)
    test_clearswi(args)
end
for args in configurations_se(["-p", phasefile_me, "-m", magfile_me, "-t", "[2,4,6]"])
    test_clearswi(args)
end

test_clearswi(["-p", phasefile_1eco, "-m", magfile_1eco, "-t", "5"])
#test_clearswi(["-p", phasefile_me_nan, "-m", magfile_me_nan_size, "-t", "[2,4]"])

## TODO: Test error and warning messages

cd(original_path)
GC.gc()
rm(tmpdir, recursive=true)
