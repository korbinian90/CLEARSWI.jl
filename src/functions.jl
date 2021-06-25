"""
    calculateSWI(data, options)

Returns the calculated SWI using 'data' and 'options'.
"""
function calculateSWI(data, options=Options())
    if !isnothing(options.writesteps) mkpath(options.writesteps) end
    getswimag(data, options) .* getswiphase(data, options)
end

createMIP(S, d=7) = createIntensityProjection(S, minimum, d)

function createIntensityProjection(S, func, d=7)
    [func(S[x,y,z:z+d-1]) for x in 1:size(S,1), y in 1:size(S,2), z in 1:size(S,3)-d+1]
end
