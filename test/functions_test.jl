using Test, ReferenceTests
using MRI, SWI

data_path = "testData/small"
TEs = [4,8,12]

mag_nii = readmag("$data_path/Mag.nii")
hdr = mag_nii.header
phase_nii = readphase("$data_path/Phase.nii")

data = Data(mag_nii, phase_nii, hdr, TEs)
o = Options()

swi = calculateSWI(data, o)

@test_reference "testData/small/swiz.sha256" swi[:,:,10]
@test_reference "testData/small/swiy.sha256" swi[:,10,:]
@test_reference "testData/small/swix.sha256" swi[10,:,:]

printstyled("Tests Passed!"; color=:green)
