data_path = "testData/small"
TEs = [4,8,12]

mag_nii = readmag("$data_path/Mag.nii")
hdr = mag_nii.header
phase_nii = readphase("$data_path/Phase.nii")

data = Data(mag_nii, phase_nii, hdr, TEs)

swi = calculateSWI(data)

options = [
    Options()
    Options(;σ=[6,6,3])
    Options(;unwrapping=:romeo)
    Options(;unwrapping=:laplacianslice)
    Options(;mode=:positive)
    Options(;mode=:negative)
    Options(;mode=:triangular)
    Options(;level=6)
    Options(;level=200)
    Options(;combination_type=:CNR=>(:wm, :gm))
    Options(;combination_type=:CNR=>(:wm, :ven))
    Options(;combination_type=:CNR=>(:wm, :csf))
    Options(;combination_type=:CNR=>(:ms, :iron_ring))
    Options(;combination_type=:SE)
    Options(;combination_type=:SE=>10)
    Options(;combination_type=:SE=>5)
]
# TODO
wrong_options = [
    Options(;level=-5)
    Options(;level=0)
    Options(;combination_type=:SE=>0)
    Options(;combination_type=:SE=>500)
]

s = []
for o in options
    push!(s, calculateSWI(data, o)[:,:,20])
end

# all results should be different
for i in 1:length(s), j in 1:(i-1)
    @test s[i] != s[j]
    if s[i] == s[j]
        @show (i,j)
    end
end
