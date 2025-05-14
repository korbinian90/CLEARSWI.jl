using QSM

# load data
data_path = CLEARSWI.dir("test", "data", "small")
TEs = [4, 8, 12]

mag_nii = readmag("$data_path/Mag.nii")
hdr = mag_nii.header
phase_nii = readphase("$data_path/Phase.nii")
data = Data(mag_nii, phase_nii, hdr, TEs)
# QSM test
qsm_run = calculateSWI(data, Options(qsm=true))
@test qsm_run != calculateSWI(data)

# QSM input test (phase in qsm)
data = Data(mag_nii, phase_nii[:,:,:,1], hdr, TEs)
qsm_run = calculateSWI(data, Options(qsm=:input))
