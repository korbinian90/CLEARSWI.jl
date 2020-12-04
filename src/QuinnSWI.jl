using DSP, ImageFiltering, MriResearchTools
function quinn(mag, phase, TEs; m=4, f=1, pth=raw"I:\Korbi\test\quinn")
    mkpath(pth)
    mag = readmag(mag)
    phase = readphase(phase)
    phase = float.(phase)
    savenii(phase[:,:,:,3], "ph3", pth)
    @time quinn_homodyne!(phase, mag, f)
    savenii(phase[:,:,:,3], "homodyne3", pth)
    @time quinn_temporal_unwrap!(phase)
    savenii(phase[:,:,:,3], "unwrapped", pth)
    @time f = MriResearchTools.calculateB0_unwrapped(phase, mag, TEs)
    savenii(f, "f", pth)
    pmask = calculate_pmask(f, TEs, :lin)
    savenii(pmask, "pmask", pth)
    avmag = mean(mag; dims=4)
    savenii(avmag, "avmag", pth)
    swi = avmag .* pmask.^m
    savenii(swi, "swi", pth)
    return swi
end

function quinn_temporal_unwrap!(phase)
    for I in CartesianIndices(phase[:,:,:,1])
        for ieco in 2:size(phase, 4)
            phase[I,ieco] = phase[I,ieco-1] + rem2pi(phase[I,ieco] - phase[I,ieco-1], RoundNearest)
        end
    end
end

function quinn_homodyne!(phase, mag, factor=1)
    for I in CartesianIndices(size(phase)[3:4])
        hann_size = round.(Int, size(phase)[1:2] .* (0.3 / factor))
        hann_window = centered(hanning(hann_size))
        slice = mag[:,:,I] .* exp.(1.0im .* phase[:,:,I])
        phase[:,:,I] = angle.(slice ./ imfilter(slice, hann_window))
    end
end

function ImageFiltering.imfilter(c::AbstractArray{<:Complex}, w)
    return complex.(imfilter(real.(c), w), imfilter(imag.(c), w))
end

function calculate_pmask(f, TEs, type=:lin)
    X = 1000 / mean(TEs) / 2 # lower X stronger masking
    L(x) = if x>X 0 elseif x<0 1 else 1 - x/X end
    H(x) = if x>X 0 elseif x<0 1 else 1/2 * (1 + cos(π*x/X)) end
    M = if type == :lin L else H end
    return float.(M.(f))
end


