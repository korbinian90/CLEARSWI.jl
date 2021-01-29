σ = [2,2,1]
unwrapping = :laplacian
mode = :linear
level = 4
tmp_folder = "tmp"
o = Options(;σ=σ, unwrapping=unwrapping, mode=mode, level=level, writedir=tmp_folder)
@test o.σ == σ
@test o.unwrapping == unwrapping
@test o.mode == mode
@test o.level == level
@test o.combination_type == :SNR
@test o.sensitivity == nothing
@test o.writedir == "tmp"
@test o.writesteps == false

mkpath(tmp_folder)
saveconfiguration(o)
text = """unwrapping: laplacian
mode: linear
level: 4
combination_type: SNR
sensitivity: nothing
writedir: tmp
writesteps: false
magscale: identity
"""
file = "$tmp_folder/settings_swi.txt"
@test read(file, String) == text
rm(file)
try rm(tmp_folder) catch end

printstyled("Tests Passed!"; color=:green)
