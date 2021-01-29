GEPCI(data::Data, options::Options) = GEPCI(data; m=options.level, writedir=options.writedir, writesteps=options.writesteps)
GEPCI(data::Data, options::NamedTuple) = GEPCI(data; options...)
function GEPCI(data::Data; m=4, writedir=raw"I:\Korbi\GEPCI\calc", overwrite=false, writesteps=false, TE_swi=20)
    TEs = data.TEs
    hdr = data.header
    mag, phase = simulate_coil_comb(data)
    
    writesteps && savenii(phase, "phase", writedir, hdr)
    writesteps && savenii(mag, "mag", writedir, hdr)

    mask = robustmask(view(mag,:,:,:,1))
    writesteps && savenii(mask, "mask", writedir, hdr)
    init(data, writedir)
    Threads.@threads for islice in 1:size(mag,3)
        if overwrite || !is_slice_finished(writedir, islice)
            println("Calculating slice: $islice")
            @time out = monoexpfit(mag, phase, mask, TEs, islice)
            save_slice(writedir, islice, out)
        end
    end
    S0, r2s, f = load_fit(writedir)

    factor = 2π * (TEs[2] .- TEs[1])
    f = romeo(f .* factor; mag=S0) ./ factor
    writesteps && savenii(f, "unwrap_f", writedir, hdr)

    f = average_filter(f)
    writesteps && savenii(f, "filtered", writedir, hdr)

    fmask = frequency_mask(f, mask)
    writesteps && savenii(fmask, "fmask", writedir, hdr)

    fmask_m = fmask.^m
    writesteps && savenii(fmask_m, "fmask_m", writedir, hdr)

    swi_like = SWI_like(S0, r2s, fmask_m, TE_swi)
    writesteps && savenii(swi_like, "swi_like", writedir, hdr)
    savenii(S0 .* exp.(-r2s .* TE_swi), "swimag", writedir, hdr)

    gepci_swi = GEPCI_SWI(r2s, fmask_m, TE_swi)
    writesteps && savenii(gepci_swi, "gepci_swi", writedir, hdr)

    gepci_t1f = GEPCI_T1F(S0, fmask_m)
    writesteps && savenii(gepci_t1f, "gepci_t1f", writedir, hdr)

    return swi_like
end

function getfns(writedir)
    return joinpath.(writedir, ["S0", "r2s", "f"] .* ".nii")
end

function init(data::Data, writedir)
    for fn in getfns(writedir)
        if !isfile(fn)
            write_emptynii(size(data.mag)[1:3], fn; header=data.header)
        end
    end
end

function init(fnmag, writedir)
    for fn in getfns(writedir)
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

function is_slice_finished(writedir, islice)
    return any(niread(first(getfns(writedir)); mmap=true, mode="r+").raw[:,:,islice] .!= 0.0)
end

function save_slice(writedir, slice, out)
    for (i, fn) in enumerate(getfns(writedir))
        niread(fn; mmap=true, mode="r+").raw[:,:,slice] .= out[i]
    end
    GC.gc()
end

function load_fit(writedir)
    return (niread(fn).raw for fn in getfns(writedir))
end

simulate_coil_comb(data::Data) = simulate_coil_comb(data.mag, data.phase)
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
            image = make_flat(mag[I,islice,:] .* exp.(1.0im .* phase[I,islice,:]))
            S0[I], r2s[I], f[I] = coef(curve_fit(model, TEs, image, p0; autodiff=:forwarddiff))
        end
    end

    ct = 0
    for I in CartesianIndices(S0)
        t = (S0[I], r2s[I], f[I])
        box = filter(ind -> checkbounds(Bool, S0, ind), getboxaround((3,3), I))
        if any(t .<= lower) || any(t .>= upper) || abs(f[I] - median(f[box])) > 0.05
            ct += 1
            
            p0 = Float64.([median(S0[box]), median(r2s[box]), median(f[box])])
            off = (0.0, 0.0, 0.01)
            lower_new = p0 .- abs.(0.05 .* p0) .- off
            upper_new = p0 .+ abs.(0.05 .* p0) .+ off
            image = make_flat(mag[I,islice,:] .* exp.(1.0im .* phase[I,islice,:]))
            S0[I], r2s[I], f[I] = coef(curve_fit(model, TEs, image, p0; lower=lower_new, upper=upper_new, autodiff=:forwarddiff))
        end
    end
    @show ct

    return S0, r2s, f
end

function monoexpfit(mag, phase, TEs)
    function model(te, p)
        res = @. p[1]^2 * exp(-p[2] * (te + TEs[1])) * exp((1.0im * 2π * p[3]) * (te - TEs[1]))
        return make_flat(res)
    end
    p0 = [0.5, 0.01, 0.01]
    lower = [0.01, 0.001, -1.0]
    upper = [5.0, 0.2, 1.0]

    image = make_flat(mag .* exp.(1.0im .* phase))
    S0, r2s, f = coef(curve_fit(model, TEs, image, p0; autodiff=:forwarddiff))

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

function frequency_mask(f, mask)
    # lower X stronger masking
    @show X = estimatequantile(f .+ (.!mask .* NaN), 0.95) / 0.95 # quantile is used because maximum is unstable
    @show X = 0.04
    L(x) = if x>X 0 elseif x<0 1 else 1 - x/X end
    return float.(L.(f))
end

function SWI_like(S0, r2s, fmask_m, TE)
    return S0 .* exp.(-r2s .* TE) .* fmask_m
end

function GEPCI_SWI(r2s, fmask_m, TE)
    return exp.(-r2s .* TE) .* fmask_m
end

function GEPCI_T1F(S0, fmask_m)
    return S0 .* fmask_m
end
