function getswiphase(data, options)
    mask = robustmask(view(data.mag,:,:,:,1))
    savenii(mask, "maskforphase", options.writesteps, data.header)
    combined = getcombinedphase(data, options, mask)
    swiphase = createphasemask!(combined, mask, options.phase_scaling_type, options.phase_scaling_strength)
    savenii(swiphase, "swiphase", options.writesteps, data.header)
    return swiphase
end

function createphasemask!(swiphase, mask, phase_scaling_type, phase_scaling_strength)
    swiphase[.!mask] .= 0
    pos = swiphase .> 0

    if phase_scaling_type == :negativetanh
        swiphase = -swiphase
        phase_scaling_type = :tanh
    end

    if phase_scaling_type == :tanh
        m = median(swiphase[mask .& (swiphase .> 0)])
        m *= 10 / phase_scaling_strength
        f(x) = (1 + tanh(1 - x/m)) / 2
        swiphase .= f.(swiphase)

    elseif phase_scaling_type == :positive
        swiphase[.!pos] .= 1
        swiphase[pos] .= robustrescale(swiphase[pos], 1, 0) .^ phase_scaling_strength

    elseif phase_scaling_type == :negative
        swiphase[pos] .= 1
        swiphase[.!pos] .= robustrescale(swiphase[.!pos], 0, 1) .^ phase_scaling_strength

    elseif phase_scaling_type == :triangular
        swiphase[pos] .= robustrescale(swiphase[pos], 1, 0) .^ phase_scaling_strength
        swiphase[.!pos] .= robustrescale(swiphase[.!pos], 0, 1) .^ phase_scaling_strength

    else
        error("$phase_scaling_type not defined for SWI")
    end
    swiphase[swiphase .< 0] .= 0
    return swiphase
end

function getcombinedphase(data, options, mask)
    save(image, name) = savenii(image, name, options.writesteps, data.header)

    if options.qsm == :input
        @assert size(data.phase, 4) == 1 "QSM input must be a single echo"
        return high_pass_qsm(data.phase, options, mask, save)

    elseif options.qsm
        return qsm_contrast(data, options, save)

    elseif options.phase_unwrap == :laplacian
        return laplacian_combine(data, options, mask, save)

    elseif options.phase_unwrap == :laplacianslice
        return laplacianslice_combine(data, options, mask, save)

    elseif options.phase_unwrap == :romeo
        return romeo_combine(data, options, mask, save)
    end

    error("Unwrapping $(options.phase_unwrap) ($(typeof(options.phase_unwrap))) not defined!")
end

function qsm_contrast(data, options, save)
    vsz = data.header.pixdim[2:4]
    TEs = to_dim(data.TEs, 4)
    mask = options.qsm_mask
    
    if isnothing(mask)
        mask = qsm_mask_filled(data.phase[:,:,:,1])
    end

    combined = qsm_romeo_B0(data.phase, data.mag, mask, TEs, vsz, B0=3; gpu=options.gpu, save, iterations=800, erosions=0)
    save(combined, "qsm_romeo_B0")
    combined = high_pass_qsm(combined, options, mask, save)
    return combined
end

function high_pass_qsm(qsm, options, mask, save)
    qsm .-= gaussiansmooth3d(qsm, options.phase_hp_sigma; mask, dims=1:2)
    save(qsm, "filteredphase")
    return qsm
end

function laplacian_combine(data, options, mask, save)
    TEs = to_dim(data.TEs, 4)
    unwrapped = similar(data.phase)
    for iEco in axes(data.phase, 4)
        unwrapped[:,:,:,iEco] .= laplacianunwrap(view(data.phase,:,:,:,iEco))
    end
    save(unwrapped, "unwrappedphase")

    for iEco in axes(data.phase, 4)
        smoothed = gaussiansmooth3d(unwrapped[:,:,:,iEco], options.phase_hp_sigma; mask, dims=1:3)
        unwrapped[:,:,:,iEco] .-= smoothed
    end
    combined = combine_phase(unwrapped, data.mag, TEs)
    save(unwrapped, "filteredphase")
    save(combined, "combinedphase")
    return combined
end

function laplacianslice_combine(data, options, mask, save)
    TEs = to_dim(data.TEs, 4)
    unwrapped = similar(data.phase)
    for iEco in axes(data.phase, 4), iSlc in axes(data.phase, 3)
        unwrapped[:,:,iSlc,iEco] .= laplacianunwrap(view(data.phase,:,:,iSlc,iEco))
        smoothed = gaussiansmooth3d(unwrapped[:,:,iSlc,iEco], options.phase_hp_sigma; mask=mask[:,:,iSlc], dims=1:2)
        unwrapped[:,:,iSlc,iEco] .-= smoothed
    end
    combined = combine_phase(unwrapped, data.mag, TEs)
    save(unwrapped, "unwrappedphase")
    save(combined, "combinedphase")
    return combined
end

function romeo_combine(data, options, mask, save)
    TEs = to_dim(data.TEs, 4)
    unwrapped = romeo(data.phase; data.mag, TEs)#, mask = mask)
    save(unwrapped, "unwrappedphase")

    combined = combine_phase(unwrapped, data.mag, TEs)
    save(combined, "combinedphase")

    combined .-= gaussiansmooth3d(combined, options.phase_hp_sigma; mask, dims=1:2)
    save(combined, "filteredphase")
    return combined
end

# returns B0 in [Hz] if TEs in [ms]
combine_phase(unwrapped::AbstractArray{T,3}, mag, TEs) where T = copy(unwrapped) # one echo
function combine_phase(unwrapped::AbstractArray{T,4}, mag, TEs) where T
    return calculateB0_unwrapped(unwrapped, mag, TEs)
end
