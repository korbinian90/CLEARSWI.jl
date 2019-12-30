using Test
using MRI, SWI

data_path = "testData/small"
TEs = [4,8,12]

mag_nii = readmag("$data_path/Mag.nii")
hdr = mag_nii.header
phase_nii = readphase("$data_path/Phase.nii")

data = Data(mag_nii, phase_nii, hdr, TEs)

swi = calculateSWI(data)

printstyled("Tests Passed!"; color=:green)
