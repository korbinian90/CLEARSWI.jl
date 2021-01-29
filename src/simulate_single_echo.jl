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
