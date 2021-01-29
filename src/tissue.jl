# https://www.sciencedirect.com/science/article/pii/S0730725X07001701?via%3Dihub
T2s = Dict( :B7T    => Dict(:gm => 33.2, :wm => 26.8, :ven => 7.4, :csf => 500, :ms => 50, :iron_ring => 21), # csf and ms estimated from 7T scan)
            :B3T    => Dict(:gm => 66.0, :wm => 53.2, :ven => 25), # ven visually estimated from swi_3T scan
            :B1p5T  => Dict(:gm => 84.0, :wm => 66.2))
#https://www.ncbi.nlm.nih.gov/pubmed/18259791
T1 = Dict(  :B7T    => Dict(:gm => 1939, :wm => 1126, :putamen => 1644, :caudate_head => 1684, #=estimated=#:ven => 2000, :ms => 1200),
            :B3T    => Dict(:gm => 1607, :wm => 838, :putamen => 1332, :caudate_head => 1395, #=estimated=#:ven => 1600, :ms => 1000),
            :B1p5T  => Dict(:gm => 1197, :wm => 646, :putamen => 1084, :caudate_head => 1109, #=estimated=#:ven => 1400, :ms => 800))
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3375320/
PD = Dict(:gm => 0.71, :wm => 0.63, :ms => 0.68, :ven => 0.9 #=estimated=#)
factor = Dict(  :B7T => Dict(:gm => 1.0, :wm => 0.95, :ven => 0.95, :csf => 0.7, :ms => 0.8, :iron_ring => 0.96), # estimated from 7T scan
                :B3T => Dict(:gm => 1.0, :wm => 0.95, :ven => 0.95)) # estimated from 7T scan

function gettissue_easy(field)
    T2s[field], factor[field]
end
