function calculateSWI(data, options)
    if !options.writesteps
        options.writedir = nothing
    end
    getswimag(data, options) .* getswiphase(data, options)
end

function getswimag(data, options)
    combined_echoes = combine_echoes_swi(data.mag, data.TEs, options.combination_type)
    @debug savenii(combined_echoes, "combined_mag", options.writedir, data.header)

    sensitivity = options.sensitivity
    if sensitivity === nothing
        sensitivity = getsensitivity(data.mag[:,:,:,1], pixdim = getpixdim(data))
    end

    @debug savenii(sensitivity, "sensitivity", options.writedir, data.header)
    combined_echoes ./ sensitivity
end

function combine_echoes_swi(mag, TEs, type)
    T2s, factor = gettissue_easy(:B7T)

    if type == :SNR
        if ndims(mag) == 3 # only one echo
            copy(mag)
        else
            RSS(mag)
        end
    elseif typeof(type) <: Pair
        type, (w1, w2) = type
        if type == :CNR
            S(TE, tissue) = factor[tissue] * exp(-TE / T2s[tissue])
            w(TE, w1, w2) = S(TE, w1) - S(TE, w2)

            combine_weighted(mag, w.(TEs, w1, w2))
        end
    else
        throw("ERROR: $type not defined for combination of echoes!")
    end
end

function combine_weighted(mag, w)
    combined = zeros(size(mag)[1:3])
    for ieco in 1:size(mag, 4)
        combined .+= mag[:,:,:,ieco] .* w[ieco]
    end
    combined
end

function getswiphase(data, options)
    mask = getrobustmask(view(data.mag,:,:,:,1))
    @debug savenii(mask, "maskforphase", options.writedir, data.header)

    combined = getcombinedphase(data, options, mask)
    scaleandthreshold!(combined, mask, options.mode, options.level)
end

function getcombinedphase(data, options, mask)
    phase = data.phase
    mag = data.mag
    TEs = data.TEs
    σ = options.σ

    # reshape TEs to the highest dimension, eg. "size(TEs) = (1, 1, 1, 3)"
    TEs = reshape(TEs, 1, 1, 1, length(TEs))
    @show typeof(phase)
    unwrapped = similar(phase)
    if options.unwrapping == :laplacian
        for iEco in 1:size(phase, 4)
            unwrapped[:,:,:,iEco] .= laplacianunwrap(view(phase,:,:,:,iEco))
        end
        @debug savenii(unwrapped, "unwrappedphase", options.writedir, data.header)

        for iEco in 1:size(phase, 4)
            unwrapped[:,:,:,iEco] .-= gaussiansmooth3d(unwrapped[:,:,:,iEco], σ; mask = mask, dims = 1:2)
        end
        combined = combine_echoes(unwrapped, mag, TEs)
        @debug savenii(unwrapped, "filteredphase", options.writedir, data.header)
        @debug savenii(combined, "combinedphase", options.writedir, data.header)

    elseif unwrapping == :laplacianslice
        for iEco in 1:size(phase, 4), iSlc in 1:size(phase, 3)
            unwrapped[:,:,iSlc,iEco] .= laplacianunwrap(view(phase,:,:,iSlc,iEco))
            unwrapped[:,:,iSlc,iEco] .-= gaussiansmooth3d(unwrapped[:,:,iSlc,iEco], σ; mask = mask[:,:,iSlc], dims = 1:2)
        end
        combined = combine_echoes(unwrapped, mag, TEs)
        @debug savenii(unwrapped, "unwrapped", options.writedir, data.header)
        @debug savenii(combined, "combinedphase", options.writedir, data.header)

    elseif unwrapping == :romeo
        unwrapped = kunwrap(phase, mag; TEs = TEs)#, mask = mask)
        @debug savenii(unwrapped, "unwrappedphase", options.writedir, data.header)

        combined = combine_echoes(unwrapped, mag, TEs)
        @debug savenii(combined, "combinedphase", options.writedir, data.header)

        combined .-= gaussiansmooth3d(combined, σ; mask = mask, dims = 1:2)
        @debug savenii(combined, "filteredphase", options.writedir, data.header)

    else
        error("Unwrapping $unwrapping not defined!")
    end

    combined
end

combine_echoes(unwrapped::AbstractArray{T,3}, mag, TEs) where T   = copy(unwrapped) # only one echo
function combine_echoes(unwrapped::AbstractArray{T,4}, mag, TEs) where T
    dim = 4
    TEs = reshape(TEs, ones(Int, dim-1)..., length(TEs)) # size = (1,1,1,nEco)

    combined = sum(unwrapped .* mag; dims = dim)
    combined ./= sum(mag .* Float32.(TEs); dims = dim)
    dropdims(combined; dims = dim)
end

function scaleandthreshold!(swiphase, mask, mode, level)
    swiphase[.!mask] .= 0
    pos = swiphase .> 0

    if mode == :tanh
        m = median(swiphase[mask .& (swiphase .> 0)])
        m *= 10 / level
        f(x) = (1 + tanh(1 - x/m)) / 2
        swiphase .= f.(swiphase)

    elseif mode == :positive
        swiphase[.!pos] .= 1
        swiphase[pos] .= robustrescale(swiphase[pos], 1, 0.5) .^ level

    elseif mode == :negative
        swiphase[pos] .= 1
        swiphase[.!pos] .= robustrescale(swiphase[.!pos], 0, 1) .^ level

    elseif mode == :triangular
        swiphase[pos] .= robustrescale(swiphase[pos], 1, 0) .^ level
        swiphase[.!pos] .= robustrescale(swiphase[.!pos], 0, 1) .^ level

    else
        error("$mode not defined for SWI")
    end
    swiphase[swiphase .< 0] .= 0
    swiphase
end

createMIP(S, d = 7) = [minimum(S[x,y,z:z+d-1]) for x in 1:size(S,1), y in 1:size(S,2), z in 1:size(S,3)-d+1]
getpixdim(data) = data.header.pixdim[2:(1+ndims(data.mag))]
