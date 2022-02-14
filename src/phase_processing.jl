function getswiphase(data, options)
    mask = robustmask(view(data.mag,:,:,:,1))
    savenii(mask, "maskforphase", options.writesteps, data.header)
    # TODO output not readable
    combined = getcombinedphase(data, options, mask)
    swiphase = creatphasemask!(combined, mask, options.phase_scaling_type, options.phase_scaling_strength)
    savenii(swiphase, "swiphase", options.writesteps, data.header)
    return swiphase
end

function creatphasemask!(swiphase, mask, phase_scaling_type, phase_scaling_strength)
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
    phase = data.phase
    mag = data.mag
    TEs = to_dim(data.TEs, 4)
    σ = options.phase_hp_σ

    unwrapped = similar(phase)
    if options.phase_unwrap == :laplacian
        for iEco in 1:size(phase, 4)
            unwrapped[:,:,:,iEco] .= laplacianunwrap(view(phase,:,:,:,iEco))
        end
        savenii(unwrapped, "unwrappedphase", options.writesteps, data.header)

        for iEco in 1:size(phase, 4)
            smoothed = gaussiansmooth3d(unwrapped[:,:,:,iEco], σ; mask, dims=1:3)
            unwrapped[:,:,:,iEco] .-= smoothed
        end
        combined = combine_phase(unwrapped, mag, TEs)
        savenii(unwrapped, "filteredphase", options.writesteps, data.header)
        savenii(combined, "combinedphase", options.writesteps, data.header)

    elseif options.phase_unwrap == :laplacianslice
        for iEco in 1:size(phase, 4), iSlc in 1:size(phase, 3)
            unwrapped[:,:,iSlc,iEco] .= laplacianunwrap(view(phase,:,:,iSlc,iEco))
            smoothed = gaussiansmooth3d(unwrapped[:,:,iSlc,iEco], σ; mask=mask[:,:,iSlc], dims=1:2)
            unwrapped[:,:,iSlc,iEco] .-= smoothed
        end
        combined = combine_phase(unwrapped, mag, TEs)
        savenii(unwrapped, "unwrappedphase", options.writesteps, data.header)
        savenii(combined, "combinedphase", options.writesteps, data.header)

    elseif options.phase_unwrap == :romeo
        unwrapped = romeo(phase; mag, TEs)#, mask = mask)
        savenii(unwrapped, "unwrappedphase", options.writesteps, data.header)

        combined = combine_phase(unwrapped, mag, TEs)
        savenii(combined, "combinedphase", options.writesteps, data.header)

        combined .-= gaussiansmooth3d(combined, σ; mask, dims=1:2)
        savenii(combined, "filteredphase", options.writesteps, data.header)

    else
        error("Unwrapping $options.phase_unwrap ($(typeof(options.phase_unwrap))) not defined!")
    end

    return combined
end

combine_phase(unwrapped::AbstractArray{T,3}, mag, TEs) where T = copy(unwrapped) # one echo
function combine_phase(unwrapped::AbstractArray{T,4}, mag, TEs) where T
    dims = 4
    TEs = reshape(TEs, ones(Int, dims-1)..., length(TEs)) # size = (1,1,1,nEco)

    combined = sum(unwrapped .* TEs .* mag .* mag; dims)
    combined ./= sum(mag .* mag .* float.(TEs).^2; dims)
    return dropdims(combined; dims)
end
