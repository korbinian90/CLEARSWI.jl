function GEPCI(fnmag, fnphase, TEs, slices=:; stored=false, TE_swi=20)
    pth = raw"I:\Korbi\GEPCI\filter"
    hdr = MriResearchTools.header(readmag(fnmag; mmap=true))
    if !stored
        mkpath(pth)
        mag = float.(readmag(fnmag))
        mask = robustmask(view(mag,:,:,:,1))
        savenii(mask, "mask", pth, hdr)
        phase = float.(readphase(fnphase))
        mag, phase = simulate_coil_comb(mag, phase)
        savenii(mag, "mag", pth, hdr)
        savenii(phase, "phase", pth, hdr)
        #temporal_unwrap!(phase)
        #romeo!(phase; mag, TEs=TEs .- TEs[1])
        #savenii(phase, "unwrapped", pth)
        init(fnmag, pth)
        Threads.@threads for islice in (1:size(mag, 3))[slices]
            if !is_slice_finished(pth, islice)
                println("Calculating slice: $islice")
                @time out = monoexpfit(mag, phase, mask, TEs, islice)
                save_slice(pth, islice, out)
            end
        end
    end
    S0, r2s, f = load_fit(pth)

    factor = 2π * (TEs[2] .- TEs[1])
    f = romeo(f .* factor; mag=S0) ./ factor
    savenii(f, "unwrap_f", pth, hdr)

    f = average_filter(f)
    savenii(f, "filtered", pth, hdr)

    fmask = frequency_mask(f)
    savenii(fmask, "fmask", pth, hdr)

    fmask4 = fmask.^4
    savenii(fmask4, "fmask4", pth, hdr)

    swi_like = SWI_like(S0, r2s, fmask, TE_swi)
    savenii(swi_like, "swi_like", pth, hdr)

    gepci_swi = GEPCI_SWI(r2s, fmask, TE_swi)
    savenii(gepci_swi, "gepci_swi", pth, hdr)

    gepci_t1f = GEPCI_T1F(S0, fmask)
    savenii(gepci_t1f, "gepci_t1f", pth, hdr)
end

function getfns(pth)
    return joinpath.(pth, ["S0", "r2s", "f"] .* ".nii")
end

function init(fnmag, pth)
    for fn in getfns(pth)
        if !isfile(fn)
            ref = niread(fnmag; mmap=true)
            write_emptynii(size(ref)[1:3], fn; header=ref.header)
        end
    end
end

function getmask(mag)
    mask = robustmask(view(mag,:,:,:,1))
    filtersize = (5,5,3)
    kernel = centered(ones(filtersize))
    kernel ./= sum(kernel)
    return imfilter(mask, kernel) .> 0.5
end

function is_slice_finished(pth, islice)
    return any(niread(first(getfns(pth)); mmap=true, mode="r+").raw[:,:,islice] .!= 0.0)
end

function save_slice(pth, slice, out)
    for (i, fn) in enumerate(getfns(pth))
        niread(fn; mmap=true, mode="r+").raw[:,:,slice] .= out[i]
    end
    GC.gc()
end

function load_fit(pth)
    return (niread(fn).raw for fn in getfns(pth))
end

function simulate_coil_comb(mag, phase)
    return mag[:,:,:,1] .* mag, rem2pi.(phase .- phase[:,:,:,1], RoundNearest)
end

function temporal_unwrap!(phase)
    for I in CartesianIndices(phase[:,:,:,1])
        for ieco in 2:size(phase, 4)
            phase[I,ieco] = phase[I,ieco-1] + rem2pi(phase[I,ieco] - phase[I,ieco-1], RoundNearest)
        end
    end
end

function monoexpfit(mag, phase, mask, TEs, islice)
    function model(te, p)
        res = @. p[1]^2 * exp(-p[2] * (te + TEs[1])) * exp((1.0im * 2π * p[3]) * (te - TEs[1]))
        return make_flat(res)
    end
    p0 = [0.5, 0.01, 0.01]
    lower = [0.01, 0.001, -1.0]
    upper = [5.0, 0.2, 1.0]

    
    S0 = zeros(Float32, size(mag)[1:2])
    r2s = zeros(Float32, size(mag)[1:2])
    f = zeros(Float32, size(mag)[1:2])

    for I in CartesianIndices(S0)
        if mask[I,islice]
            image = make_flat(view(mag,I,islice,:) .* exp.(1.0im .* view(phase,I,islice,:)))
            S0[I], r2s[I], f[I] = coef(curve_fit(model, TEs, image, p0; lower, upper, autodiff=:forwarddiff))
        end
    end

    ct = 0
    for I in CartesianIndices(S0)
        t = (S0[I], r2s[I], f[I])
        box = filter(ind -> checkbounds(Bool, S0, ind) && ind != I, getboxaround((3,3), I))
        if any(t .== lower) || any(t .== upper) || abs(f[I] - median(f[box])) > 0.1
            ct += 1
            
            p0 = Float64.([median(S0[box]), median(r2s[box]), median(f[box])])
            off = (0.0, 0.0, 0.01)
            lower_new = p0 .- abs.(0.05 .* p0) .- off
            upper_new = p0 .+ abs.(0.05 .* p0) .+ off
            image = make_flat(view(mag,I,islice,:) .* exp.(1.0im .* view(phase,I,islice,:)))
            S0[I], r2s[I], f[I] = coef(curve_fit(model, TEs, image, p0; lower=lower_new, upper=upper_new, autodiff=:forwarddiff))
        end
    end
    @show ct

    return S0, r2s, f
end

function getboxaround(boxsize, idx=CartesianIndex(zeros(Int, length(boxsize))...))
    width = div.((boxsize .- 1), 2)
    r(x) = -x:x
    return idx .+ CartesianIndices(r.(width))
end

function make_flat(c)
    return [real.(c)..., imag.(c)...]
end

function average_filter(image; size=(7,7,1))
    kernel = centered(ones(size))
    kernel ./= sum(kernel) 
    return image .- imfilter(image, kernel)
end

function frequency_mask(f)
    @show X = maximum(f) # lower X stronger masking
    @show X = estimatequantile(f, 0.95)
    L(x) = if x>X 0 elseif x<0 1 else 1 - x/X end
    return float.(L.(f))
end

function SWI_like(S0, r2s, fmask, TE)
    return S0 .* exp.(-r2s .* TE) .* fmask.^4
end

function GEPCI_SWI(r2s, fmask, TE)
    return exp.(-r2s .* TE) .* fmask.^4
end

function GEPCI_T1F(S0, fmask)
    return S0 .* fmask.^4
end
