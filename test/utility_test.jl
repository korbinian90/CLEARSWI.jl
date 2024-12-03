σ = [2,2,1]
phase_unwrap = :laplacian
phase_scaling_type = :linear
phase_scaling_strength = 4
tmp_folder = "tmp"
mag_combine = (:CNR=>(:gm, :wm))
o = Options(mag_combine=mag_combine, phase_hp_sigma=σ, phase_unwrap=phase_unwrap, phase_scaling_type=phase_scaling_type, phase_scaling_strength=phase_scaling_strength, writesteps=tmp_folder)
@test o.phase_hp_sigma == σ
@test o.phase_unwrap == phase_unwrap
@test o.phase_scaling_type == phase_scaling_type
@test o.phase_scaling_strength == phase_scaling_strength
@test o.mag_combine == (:CNR=>(:gm, :wm))
@test o.mag_sens === nothing
@test o.writesteps == "tmp"
@test o.mag_softplus == true

mkpath(tmp_folder)
saveconfiguration(o)
text = """mag_combine: :CNR => (:gm, :wm)
mag_sens: nothing
mag_softplus: true
phase_unwrap: laplacian
phase_hp_sigma: [2, 2, 1]
phase_scaling_type: linear
phase_scaling_strength: 4
writesteps: tmp
qsm: false
CLEARSWI.jl github version-tag: $(pkgversion(CLEARSWI))
"""
file = "$tmp_folder/settings_swi.txt"
@test read(file, String) == text
rm(file)
try rm(tmp_folder) catch end

@test isdir(CLEARSWI.dir("test", "data", "small"))
