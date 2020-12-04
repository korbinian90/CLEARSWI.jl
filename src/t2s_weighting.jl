function simulate_late_mag(mag, TEs, sim_te)
    @time t2s = NumART2star(mag, TEs)
    t2s[t2s .> 100] .= 100
    t2s[t2s .< 1] .= 1
    t2s[.!isfinite.(t2s)] .= 100
    #@time gaussiansmooth3d!(t2s; nbox=2, boxsizes=[[3,3], [3,3], [1,1]])
    @time sens = getsensitivity(mag)
    # TODO add small average filter (gaussian 2box 3x3x1)
    #? sim_te = sim_te - weighted_mean_te?
    return RSS(mag) ./ sens .* exp.(-sim_te ./ t2s)
end
