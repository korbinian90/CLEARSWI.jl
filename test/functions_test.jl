# load data
data_path = CLEARSWI.dir("test", "data", "small")
TEs = [4, 8, 12]

mag_nii = readmag("$data_path/Mag.nii")
hdr = mag_nii.header
phase_nii = readphase("$data_path/Phase.nii")

m = "Multi-echo data needs the echo times (TEs)"
@test_throws ErrorException(m) Data(mag_nii, phase_nii, hdr)

data = Data(mag_nii, phase_nii, hdr, TEs)

# echo formatting test for 2d row vector
Data(mag_nii, phase_nii, hdr, [4 8 12])

# default test
swi = calculateSWI(data, Options(phase_unwrap=:romeo))
mip = createMIP(swi)
mIP = createIntensityProjection(swi, minimum)
MIP = createIntensityProjection(swi, maximum)
meanIP = createIntensityProjection(swi, mean)

# single-echo
se_data = Data(mag_nii[:,:,:,1], phase_nii[:,:,:,1], hdr)
for o in (Options(), Options(phase_unwrap=:romeo))
    se_swi = calculateSWI(se_data, o)
    @test size(se_swi) == size(mag_nii)[1:3]
    @test all(isfinite, se_swi)
    createMIP(se_swi)
end

options = [
    Options()
    Options(phase_hp_sigma=[6,6,3])
    Options(phase_scaling_type=:positive)
    Options(phase_scaling_type=:negativetanh)
    Options(phase_scaling_type=:negative)
    Options(phase_scaling_type=:triangular)
    Options(phase_scaling_strength=6)
    Options(phase_scaling_strength=200)
    Options(mag_combine=:CNR=>(:wm, :gm))
    Options(mag_combine=:CNR=>(:wm, :ven))
    Options(mag_combine=:CNR=>(:wm, :csf))
    Options(mag_combine=:CNR=>(:ms, :iron_ring))
    Options(mag_combine=:average)
    Options(mag_combine=:SE=>10)
    Options(mag_combine=:SE=>5)
    Options(mag_combine=(:closest=>5)) # 1
    Options(mag_combine=(:echo=>2)) # 2
    Options(mag_combine=:last) # 3
    Options(mag_softplus=false)
    Options(mag_sens=[1])
    Options(mag_sens=(:sigma_mm => 3))
]
# TODO
wrong_options = [
    Options(;phase_scaling_strength=-5)
    Options(;phase_scaling_strength=0)
    Options(;mag_combine=:SE=>0)
    Options(;mag_combine=:SE=>500)
]

s = []
for o in options
    # change phase unwrap to romeo
    kw = Dict(f=>getfield(o,f) for f in fieldnames(typeof(o)))
    o = Options(; kw..., phase_unwrap=:romeo)
    push!(s, calculateSWI(data, o)[:,:,20])
end

# all results should be different
for i in eachindex(s), j in 1:(i-1)
    @test s[i] != s[j]
    if s[i] == s[j]
        @show (i,j)
    end
end

# a mag_sens pair other than :sigma_mm => value is an error, not a division by the pair
@test_throws ArgumentError calculateSWI(data, Options(mag_sens=:sigma=>7))

# Without a QSM backend, which a statically compiled program leaves out, calculateSWI
# infers as concrete arrays. A loaded backend adds the QSM path, whose type is not
# inferable, so the check depends on no other test item having loaded one.
qsm_backend = any(ext -> Base.get_extension(CLEARSWI.MriResearchTools, ext) !== nothing,
                  (:QSMExt, :QuantitativeSusceptibilityMappingTGVExt))
for o in (Options(), Options(mag_sens=:sigma_mm=>7), Options(phase_scaling_type=:positive, phase_scaling_strength=2))
    T = only(Base.return_types(calculateSWI, (typeof(data), typeof(o))))
    @test all(isconcretetype, Base.uniontypes(T)) skip=qsm_backend
end
