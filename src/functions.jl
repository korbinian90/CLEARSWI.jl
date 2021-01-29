
"""
    calculateSWI(data, options)

Returns the calculated SWI using 'data' and 'options'.
"""
function calculateSWI(data, options=Options())
    if !options.writesteps
        options.writedir = nothing
    end
    getswimag(data, options) .* getswiphase(data, options)
end

function getswimag(data, options)
    combined_echoes = combine_echoes_swi(data.mag, data.TEs, options.combination_type)
    savenii(combined_echoes, "combined_mag", options.writedir, data.header)

    sensitivity = options.sensitivity
    if sensitivity === nothing
        sensitivity = getsensitivity(data.mag[:,:,:,1], getpixdim(data))
    end
    savenii(sensitivity, "sensitivity", options.writedir, data.header)

    swimag = options.magscale(combined_echoes ./ sensitivity)
    savenii(swimag, "swimag", options.writedir, data.header)
    return swimag
end

function softplus(val, offset, factor=1)
    f = factor / offset
    sp(x) = log(1 + exp(f * (x - offset))) / f
    return sp(val) - sp(0)
end

function combine_echoes_swi(mag, TEs, type)
    if ndims(mag) == 3 # only one echo
        return copy(mag)
    elseif type == :SNR
        return RSS(mag)
    elseif type == :SE
        return simulate_single_echo_mag(mag, TEs)
    elseif type == :average
        return sum(mag; dims=4)
    elseif type == :last
        return mag[:,:,:,end]
    elseif typeof(type) <: Pair
        type, para = type
        if type == :CNR
            (w1, w2) = para
            weighting = calculate_cnr_weighting(TEs, w1, w2)
            return combine_weighted(mag, weighting)
        elseif type == :SE
            TE_SE = para
            return simulate_single_echo_mag(mag, TEs, TE_SE)
        elseif type == :closest
            TE_SE = para
            eco = findmin(abs.(TEs .- TE_SE))[2]
            return mag[:,:,:,eco]
        elseif type == :echo
            return mag[:,:,:,para]
        end
    end
    throw("ERROR: $type not defined for combination of echoes!")
end

function calculate_cnr_weighting(TEs, w1, w2; field=:B7T)
    T2s, factor = gettissue_easy(field)
    S(TE, tissue) = factor[tissue] * exp(-TE / T2s[tissue])
    w(TE, w1, w2) = S(TE, w1) - S(TE, w2)
    return w.(TEs, w1, w2)
end

function combine_weighted(mag, w)
    combined = zeros(size(mag)[1:3])
    for ieco in 1:size(mag, 4)
        combined .+= mag[:,:,:,ieco] .* w[ieco]
    end
    combined
end

function getswiphase(data, options)
    mask = robustmask(view(data.mag,:,:,:,1))
    savenii(mask, "maskforphase", options.writedir, data.header)
    # TODO output not readable
    combined = getcombinedphase(data, options, mask)
    swiphase = scaleandthreshold!(combined, mask, options.mode, options.level)
    savenii(swiphase, "swiphase", options.writedir, data.header)
    return swiphase
end

function getcombinedphase(data, options, mask)
    phase = data.phase
    mag = data.mag
    TEs = to_dim(data.TEs, 4)
    σ = options.σ

    unwrapped = similar(phase)
    if options.unwrapping == :laplacian
        for iEco in 1:size(phase, 4)
            unwrapped[:,:,:,iEco] .= laplacianunwrap(view(phase,:,:,:,iEco))
        end
        savenii(unwrapped, "unwrappedphase", options.writedir, data.header)

        for iEco in 1:size(phase, 4)
            smoothed = gaussiansmooth3d(unwrapped[:,:,:,iEco], σ; mask=mask, dims=1:3)
            unwrapped[:,:,:,iEco] .-= smoothed
        end
        combined = combine_echoes(unwrapped, mag, TEs)
        savenii(unwrapped, "filteredphase", options.writedir, data.header)
        savenii(combined, "combinedphase", options.writedir, data.header)

    elseif options.unwrapping == :laplacianslice
        for iEco in 1:size(phase, 4), iSlc in 1:size(phase, 3)
            unwrapped[:,:,iSlc,iEco] .= laplacianunwrap(view(phase,:,:,iSlc,iEco))
            smoothed = gaussiansmooth3d(unwrapped[:,:,iSlc,iEco], σ; mask=mask[:,:,iSlc], dims=1:2)
            unwrapped[:,:,iSlc,iEco] .-= smoothed
        end
        combined = combine_echoes(unwrapped, mag, TEs)
        savenii(unwrapped, "unwrapped", options.writedir, data.header)
        savenii(combined, "combinedphase", options.writedir, data.header)

    elseif options.unwrapping == :romeo
        unwrapped = romeo(phase, mag=mag, TEs=TEs)#, mask = mask)
        savenii(unwrapped, "unwrappedphase", options.writedir, data.header)

        combined = combine_echoes(unwrapped, mag, TEs)
        savenii(combined, "combinedphase", options.writedir, data.header)

        combined .-= gaussiansmooth3d(combined, σ; mask = mask, dims = 1:2)
        savenii(combined, "filteredphase", options.writedir, data.header)

    else
        error("Unwrapping $unwrapping ($(typeof(unwrapping))) not defined!")
    end

    combined
end

combine_echoes(unwrapped::AbstractArray{T,3}, mag, TEs) where T = copy(unwrapped) # one echo
function combine_echoes(unwrapped::AbstractArray{T,4}, mag, TEs) where T
    dims = 4
    TEs = reshape(TEs, ones(Int, dims-1)..., length(TEs)) # size = (1,1,1,nEco)

    combined = sum(unwrapped .* TEs .* mag .* mag; dims=dims)
    combined ./= sum(mag .* mag .* float.(TEs).^2; dims=dims)
    dropdims(combined; dims=dims)
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

function createMIP(S, d=7)
    [minimum(S[x,y,z:z+d-1]) for x in 1:size(S,1), y in 1:size(S,2), z in 1:size(S,3)-d+1]
end
getpixdim(data) = data.header.pixdim[2:(1+ndims(data.mag))]

function simulate_single_echo_mag(mag, TEs, TE_SE=mean(TEs))
    weighting = get_single_echo_weighting(TEs, TE_SE)
    weighting = to_dim(weighting, 4)
    return dropdims(sum(mag .* weighting; dims=4); dims=4)
end

function get_single_echo_weighting(TEs, TE_SE)
    if TE_SE < TEs[1] || TE_SE > TEs[end]
        error("Not possible to simulate TE=$TE_SE from $TEs !")
    end
    ΔTE = TEs[2] - TEs[1]
    echoend_ME = TEs[end] + ΔTE / 2
    echostart_ME = TEs[1] - ΔTE / 2
    echowidth_sim = 2min(echoend_ME - TE_SE, TE_SE - echostart_ME)

    echostart_sim = TE_SE - echowidth_sim / 2
    echoend_sim = TE_SE + echowidth_sim / 2

    return get_single_echo_weighting(TEs, echostart_sim, echoend_sim)
end

function get_single_echo_weighting(TEs, echostart_sim, echoend_sim)
    ΔTE = TEs[2] - TEs[1]
    echoend_ME = TEs[end] + ΔTE / 2
    echostart_ME = TEs[1] - ΔTE / 2
    if (echostart_sim + 1e-5) < echostart_ME || (echoend_sim - 1e-5) > echoend_ME || (echoend_sim - echostart_sim + 1e-5) < ΔTE
        error("Not possible to simulate [$(echostart_sim);$(echoend_sim)] from $TEs !")
    end
    weighting = ones(length(TEs))
    if echostart_sim > echostart_ME
        lowechoborder = argmin(abs.(TEs .- echostart_sim))
        weighting[1:lowechoborder] .= 0
        weighting[lowechoborder] = (TEs[lowechoborder] + ΔTE/2 - echostart_sim) / ΔTE
    elseif echoend_sim < echoend_ME
        highechoborder = argmin(abs.(TEs .- echoend_sim))
        weighting[highechoborder:end] .= 0
        weighting[highechoborder] = (echoend_sim - (TEs[highechoborder] - ΔTE/2)) / ΔTE
    end
    return weighting
end
