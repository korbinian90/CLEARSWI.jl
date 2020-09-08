struct Data
    mag::AbstractArray
    phase::AbstractArray
    header
    TEs::AbstractArray
end

mutable struct Options
    σ::AbstractArray
    unwrapping::Symbol
    mode::Symbol
    level::Number
    combination_type::Union{Symbol, Pair{Symbol,<:Number}, Pair{Symbol,Tuple{Symbol,Symbol}}}
    sensitivity::Union{AbstractArray, Nothing}
    writedir::Union{AbstractString, Nothing}
    writesteps::Bool
end

function Options(;σ=[4,4,0], unwrapping=:laplacian, mode=:tanh, level=4, combination_type=:SNR, sensitivity=nothing, writedir=nothing, writesteps=false)
    Options(σ, unwrapping, mode, level, combination_type, sensitivity, writedir, writesteps)
end

Base.string(t::Pair{Symbol,Tuple{Symbol,Symbol}}) = "$(t[1])_$(t[2][1])_$(t[2][2])"
Base.string(s::Nothing) = "nothing"
function saveconfiguration(options)
    open(joinpath(options.writedir, "settings_swi.txt"), "w") do io
        for fname in fieldnames(typeof(options))
            val = getfield(options, fname)
            if !(typeof(val) <: AbstractArray)
                println(io, "$fname: " * string(val))
            end
        end
    end
end
