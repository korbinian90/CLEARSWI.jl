# phase_unwrap_kernel selects the Laplacian discretisation behind
# phase_unwrap=:laplacian / :laplacianslice:
#   :dct (default) -> MriResearchTools.laplacianunwrap   (Schofield-Zhu, mirror boundaries)
#   :fft           -> MriResearchTools.laplacianunwrap_fft (discrete k-space stencil, periodic; Bilgic)
# This pins the default, checks both kernels run end-to-end on real data, and
# asserts the two discretisations produce different SWI (i.e. the knob is not a
# no-op). It also exercises the slicewise (2D) path and rejects unknown kernels.

data_path = CLEARSWI.dir("test", "data", "small")
TEs = [4, 8, 12]
mag_nii = readmag("$data_path/Mag.nii")
hdr = mag_nii.header
phase_nii = readphase("$data_path/Phase.nii")
data = Data(mag_nii, phase_nii, hdr, TEs)

# default kernel is :dct, and the field round-trips through the constructor
@test Options().phase_unwrap_kernel == :dct
@test Options(phase_unwrap_kernel=:fft).phase_unwrap_kernel == :fft

# both kernels run through the full 3D :laplacian pipeline and give finite SWI
swi_dct = calculateSWI(data, Options(phase_unwrap=:laplacian, phase_unwrap_kernel=:dct))
swi_fft = calculateSWI(data, Options(phase_unwrap=:laplacian, phase_unwrap_kernel=:fft))
@test size(swi_dct) == size(swi_fft)
@test all(isfinite, swi_dct)
@test all(isfinite, swi_fft)

# the two discretisations must actually differ (otherwise the knob does nothing)
@test swi_dct != swi_fft

# the slicewise (2D) unwrap path also honours the kernel selection
swi_slice_fft = calculateSWI(data, Options(phase_unwrap=:laplacianslice, phase_unwrap_kernel=:fft))
@test all(isfinite, swi_slice_fft)

# an unknown kernel is rejected with a clear error
@test_throws ErrorException calculateSWI(data, Options(phase_unwrap=:laplacian, phase_unwrap_kernel=:wavelet))
