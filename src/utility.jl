struct Data
    mag::AbstractArray
    phase::AbstractArray
    header
    TEs::AbstractVector
    Data(mag, phase, header, TEs=1:size(mag,4)) = new(mag, phase, header, vec(TEs))
end

struct Options
    mag_combine
    mag_sens::Union{AbstractArray, Nothing}
    mag_softplus
    phase_unwrap::Symbol
    phase_hp_σ::AbstractArray
    phase_scaling_type::Symbol
    phase_scaling_strength::Real
    writesteps::Union{AbstractString, Nothing}
end

function Options(; mag_combine=:SNR, mag_sens=nothing, mag_softplus=true, phase_unwrap=:laplacian, phase_hp_σ=[4,4,0], phase_scaling_type=:tanh, phase_scaling_strength=4, writesteps=nothing)
    Options(mag_combine, mag_sens, mag_softplus, phase_unwrap, phase_hp_σ, phase_scaling_type, phase_scaling_strength, writesteps)
end

function saveconfiguration(options)
    open(joinpath(options.writesteps, "settings_swi.txt"), "w") do io
        for fname in fieldnames(typeof(options))
            val = getfield(options, fname)
            if !(val isa AbstractArray && !(val isa Vector))
                println(io, "$fname: " * string(val))
            end
        end
    end
end
