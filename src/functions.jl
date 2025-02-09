"""
    calculateSWI(data::Data, options::Options=Options())

Returns the calculated SWI using 'data' and 'options'.

# Examples
```julia-repl
julia> TEs = [4,8,12]
julia> data = Data(mag, phase, header(mag), TEs);
julia> swi = calculateSWI(data);
```
# With Options
$(@doc Options)

# Examples
```julia-repl
julia> TEs = [4,8,12]
julia> data = Data(mag, phase, header(mag), TEs);
julia> options = Options(phase_hp_sigma=[10,10,5], mag_softplus=false)
julia> swi = calculateSWI(data, options);
```
"""
function calculateSWI(data::Data, options::Options=Options())
    getswimag(data, options) .* getswiphase(data, options)
end

"""
    createMIP(S::AbstractArray{<:Number,3}, d=7)

Creates minimum intensity projection of `S` over `d` slices.

# Examples
```julia-repl
julia> TEs = [4,8,12]
julia> data = Data(mag, phase, header(mag), TEs);
julia> swi = calculateSWI(data);
julia> mip = createMIP(swi);
```

See also [`createIntensityProjection`](@ref)
"""
createMIP(S::AbstractArray{<:Number,3}, d=7) = createIntensityProjection(S, minimum, d)


"""
    createIntensityProjection(S::AbstractArray{T,3}, func, d=7) where T

Calculates an intensity projection over the 3rd dimension of `S` with the function `func` over `d` slices.
Good function are `minimum`, `maximum`, `mean`, `std`, ...

# Examples
```julia-repl
julia> using Statistics
julia> a = rand(Float64, (64,64,20))
julia> createIntensityProjection(a, std)
```
"""
function createIntensityProjection(S::AbstractArray{T,3}, func, d=7) where T
    [func(S[x,y,z:z+d-1]) for x in axes(S,1), y in axes(S,2), z in 1:size(S,3)-d+1]
end
