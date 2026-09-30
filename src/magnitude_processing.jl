# The options are passed on one by one rather than as `Options`, so that each
# step is compiled once per type of the options it uses, not once per
# combination of all of them.
getswimag(data, options) = getswimag(data, options.mag_combine, options.mag_sens, options.mag_softplus, options.writesteps)

function getswimag(data, mag_combine, mag_sens, mag_softplus, writesteps)
    combined_mag = combine_echoes_swi(data.mag, data.TEs, mag_combine)
    savenii(combined_mag, "combined_mag", writesteps, data.header)
    swimag = sensitivity_correction(combined_mag, data, mag_sens, writesteps)
    if mag_softplus !== false
        savenii(swimag, "sensitivity_corrected_mag", writesteps, data.header)
        swimag = softplus_scaling(swimag, mag_softplus)
    end
    savenii(swimag, "swimag", writesteps, data.header)
    return swimag
end

function sensitivity_correction(combined_mag, data, sensitivity, writesteps)
    if isnothing(sensitivity)
        sensitivity = getsensitivity(data.mag, getpixdim(data))
    elseif sensitivity isa Pair
        first(sensitivity) == :sigma_mm || throw(ArgumentError("mag_sens must be an array, nothing or :sigma_mm => value"))
        sensitivity = getsensitivity(data.mag, getpixdim(data); sigma_mm=last(sensitivity))
    end
    savenii(sensitivity, "sensitivity", writesteps, data.header)
    return combined_mag ./ sensitivity
end
getpixdim(data) = ntuple(i -> data.header.pixdim[i+1], Val(ndims(data.mag)))

combine_echoes_swi(mag, TEs, type) = ndims(mag) == 3 ? copy(mag) : combine_echoes(mag, TEs, type) # 3D is only one echo

function combine_echoes(mag, TEs, type::Symbol)
    type == :SNR && return RSS(mag)
    type == :average && return dropdims(sum(mag; dims=4); dims=4)
    type == :last && return mag[:,:,:,end]
    error("$type not defined for combination of echoes!")
end
function combine_echoes(mag, TEs, type::Pair)
    kind, para = type
    kind == :CNR && para isa Union{Tuple,AbstractVector} && return combine_cnr(mag, TEs, para)
    kind == :SE && return simulate_single_echo_mag(mag, TEs, para)
    kind == :closest && return mag[:,:,:,findmin(abs.(TEs .- para))[2]]
    kind == :echo && return mag[:,:,:,echo_index(para)]
    kind == :average && return sum(mag[:,:,:,echo_index(para)]; dims=4)
    error("$kind not defined for combination of echoes!")
end
combine_echoes(mag, TEs, type) = throw(ArgumentError("mag_combine must be a Symbol or a Pair"))

# The parameter of an echo combination that selects echoes. A number that is
# not an integer is an error rather than an index, which is also what lets a
# compiled program leave out indexing with a float.
echo_index(i::Integer) = i
echo_index(x::Real) = isinteger(x) ? Int(x) : throw(ArgumentError("echo $x is not an integer"))
echo_index(v::AbstractVector) = v

function combine_cnr(mag, TEs, para)
    if length(para) == 2
        (w1, w2) = para
        field = :B7T
    else
        (w1, w2, field) = para
    end
    weighting = calculate_cnr_weighting(TEs, w1, w2; field=field)
    return combine_weighted(mag, weighting)
end

function calculate_cnr_weighting(TEs, w1::Tuple, w2::Tuple; unusedargs...)
    S(TE, w) = w[1] * exp(-TE / w[2])
    w(TE, w1, w2) = S(TE, w1) - S(TE, w2)
    return w.(TEs, Ref(w1), Ref(w2))
end

function calculate_cnr_weighting(TEs, w1, w2; field=:B7T)
    T2s, factor = gettissue_easy(field)
    S(TE, tissue) = factor[tissue] * exp(-TE / T2s[tissue])
    w(TE, w1, w2) = S(TE, w1) - S(TE, w2)
    return w.(TEs, w1, w2)
end

function combine_weighted(mag, w)
    combined = zeros(size(mag)[1:3])
    for ieco in axes(mag, 4)
        combined .+= mag[:,:,:,ieco] .* w[ieco]
    end
    return combined ./ sum(w)
end

function softplus_scaling(mag, para)
    q = estimatequantile(mag, 0.8)
    if para === true
        return softplus.(mag, q/2)
    elseif para isa Number
        return softplus.(mag, para * q)
    elseif para isa Tuple
        return softplus.(mag, para[1] * q, para[2])
    else
        error("wrong input for softplus: got $para")
    end
end

function softplus(val, offset, factor=2)
    f = factor / offset
    # stable implementation
    function sp(x)
        arg = f * (x - offset)
        soft = log(1 + exp(-abs(arg))) + max(0, arg)
        return soft / f
    end
    return sp(val) - sp(0)
end
