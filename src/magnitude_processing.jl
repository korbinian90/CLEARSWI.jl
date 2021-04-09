function getswimag(data, options)
    combined_mag = combine_echoes_swi(data.mag, data.TEs, options.mag_combine)
    savenii(combined_mag, "combined_mag", options.writesteps, data.header)
    swimag = sensitivity_correction(combined_mag, data, options)
    if options.mag_softplus
        savenii(swimag, "sensitivity_corrected_mag", options.writesteps, data.header)
        swimag = softplus_scaling(swimag)
    end
    savenii(swimag, "swimag", options.writesteps, data.header)
    return swimag
end

function sensitivity_correction(combined_mag, data, options)
    sensitivity = options.mag_sens
    if isnothing(sensitivity)
        sensitivity = getsensitivity(data.mag[:,:,:,1], getpixdim(data))
    end
    savenii(sensitivity, "sensitivity", options.writesteps, data.header)
    return combined_mag ./ sensitivity
end
getpixdim(data) = data.header.pixdim[2:(1+ndims(data.mag))]

function combine_echoes_swi(mag, TEs, type)
    if ndims(mag) == 3 # only one echo
        return copy(mag)
    elseif type == :SNR
        return RSS(mag)
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
    return combined ./ sum(w)
end

function softplus_scaling(mag)
    q = estimatequantile(mag, 0.8)
    softplus.(mag, q/2)
end

function softplus(val, offset, factor=2)
    f = factor / offset
    sp(x) = log(1 + exp(f * (x - offset))) / f
    return sp(val) - sp(0)
end
