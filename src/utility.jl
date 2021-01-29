struct Data
    mag::AbstractArray
    phase::AbstractArray
    header
    TEs::AbstractArray
end

mutable struct Options
    mag_combine::Union{Symbol, Pair{Symbol,<:Real}, Pair{Symbol,Tuple{Symbol,Symbol}}}
    mag_sens::Union{AbstractArray, Nothing}
    mag_softplus::Bool
    phase_unwrap::Symbol
    phase_hp_σ::AbstractArray
    phase_scaling_type::Symbol
    phase_scaling_strength::Real
    writesteps::Union{AbstractString, Nothing}
end

function Options(; mag_combine=:SNR, mag_sens=nothing, mag_softplus=true, phase_unwrap=:laplacian, phase_hp_σ=[4,4,0], phase_scaling_type=:tanh, phase_scaling_strength=4, writesteps=nothing)
    Options(mag_combine, mag_sens, mag_softplus, phase_unwrap, phase_hp_σ, phase_scaling_type, phase_scaling_strength, writesteps)
end

makestring(t::Pair{Symbol,Tuple{Symbol,Symbol}}) = "$(t[1])_$(t[2][1])_$(t[2][2])"
makestring(s::Nothing) = "nothing"
makestring(s) = string(s)
function saveconfiguration(options)
    open(joinpath(options.writesteps, "settings_swi.txt"), "w") do io
        for fname in fieldnames(typeof(options))
            val = getfield(options, fname)
            if !(typeof(val) <: AbstractArray)
                println(io, "$fname: " * makestring(val))
            end
        end
    end
end
