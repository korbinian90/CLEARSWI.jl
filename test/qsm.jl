@testitem "qsm" begin
cd(@__DIR__)
# load data
data_path = "data/small"
TEs = [4, 8, 12]

mag_nii = readmag("$data_path/Mag.nii")
hdr = mag_nii.header
phase_nii = readphase("$data_path/Phase.nii")
data = Data(mag_nii, phase_nii, hdr, TEs)
# QSM test
qsm_run = calculateSWI(data, Options(qsm=true))
@test qsm_run != calculateSWI(data)
end
