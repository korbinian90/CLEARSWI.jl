using QuantitativeSusceptibilityMappingTGV

# load data
data_path = CLEARSWI.dir("test", "data", "small")
TEs = [4, 8, 12]

mag_nii = readmag("$data_path/Mag.nii")
hdr = mag_nii.header
phase_nii = readphase("$data_path/Phase.nii")
data = Data(mag_nii, phase_nii, hdr, TEs)
# The QSM path. CLEARSWI reaches it through MriResearchTools.qsm_B0, which is
# answered by whichever backend extension is loaded; TGV is the one the shipped
# app depends on, so it is the one tested here. There used to be a second,
# identical test item that imported QSM.jl instead - see the note in test/Project.toml.
qsm_run = calculateSWI(data, Options(qsm=true))
@test qsm_run != calculateSWI(data)
@test size(qsm_run) == size(mag_nii)[1:3]

# 50 of 32000 voxels (0.16%) come out NaN on this dataset. Recorded, not fixed:
# `high_pass_qsm` subtracts `gaussiansmooth3d(qsm, [4,4,0]; mask, dims=1:2)` from
# the whole array, and the masked smoother returns NaN wherever an in-plane
# window contains no masked voxel at all (0/0). All 50 sit *outside* the QSM
# mask, but they still reach `filteredphase`, `swiphase` and the SWI, and they
# are not background air - the magnitude there averages 45% of the image maximum,
# because the phase-quality QSM mask excludes voxels that carry real signal.
# Fixing it means either confining the subtraction to the mask or giving the
# masked smoother a no-support fallback, both of which change numbers in a
# shipped path, so it is the maintainer's call.
@test_broken all(isfinite, qsm_run)

# a pre-calculated QSM supplied as the phase input
data_input = Data(mag_nii, phase_nii[:,:,:,1], hdr, TEs)
qsm_input_run = calculateSWI(data_input, Options(qsm=:input))
@test size(qsm_input_run) == size(mag_nii)[1:3]
@test all(isfinite, qsm_input_run)
