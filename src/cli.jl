# The clearswi command line. It parses with MriResearchTools.CLI and reads NIfTI
# through loaders of fixed type, so that it compiles with `juliac --trim` as well
# as running in Julia.

const CLI = MriResearchTools.CLI

const OPTIONS = [
    CLI.Option("--magnitude", "-m", "The magnitude image (single or multi-echo)"),
    CLI.Option("--phase", "-p", "The phase image (single or multi-echo)"),
    CLI.Option("--output", "-o", "The output path or filename (default: clearswi.nii)"),
    CLI.Option("--echo-times", "-t", """The echo times are required for multi-echo datasets
        specified in array or range syntax (eg. "[1.5,3.0]" or
        "3.5:3.5:14")."""; nargs=:many),
    CLI.Option("--mip-slices", "-s", "The number of slices in the MIP image (default: 7)"),
    CLI.Option("--qsm", "", "When activated uses TGV QSM for phase weighting."; nargs=:none),
    CLI.Option("--qsm-input", "", """Give pre-calculated QSM instead of phase as input. Not fine-tuned!
        Please adjust phase filter to something like:
        "--filter-size [20,20,0] --phase-phase_scaling_strength 0.04" """),
    CLI.Option("--qsm-mask", "", """The mask used for QSM. Use a custom mask, if the qsm_mask.nii
        is not good for your data."""),
    CLI.Option("--mag-combine", "", """SNR | average | echo <n> | SE <te>.
        Magnitude combination algorithm. echo <n> selects a specific
        echo; SE <te> simulates a single echo scan of the given echo
        time. (default: SNR)"""; nargs=:many),
    CLI.Option("--mag-sensitivity-correction", "", """ <filename> | on | off.
        Use the CLEAR-SWI sensitivity correction. Alternatively, a
        sensitivity map can be read from a file (default: on)"""),
    CLI.Option("--mag-softplus-scaling", "", """on | off.
        Set softplus scaling of the magnitude (default: on)"""),
    CLI.Option("--unwrapping-algorithm", "", "laplacian | romeo | laplacianslice (default: laplacian)"),
    CLI.Option("--filter-size", "", """Size for the high-pass phase filter in voxels. Can be
        given as <x> <y> <z> or in array syntax (e.g. [2.2,3.1,0],
        which is effectively a 2D filter). (default: [4,4,0])"""; nargs=:many),
    CLI.Option("--phase-scaling-type", "", """tanh | negativetanh | positive | negative | triangular
        Select the type of phase scaling. positive or negative with a
        strength of 3-6 is used in standard SWI. (default: tanh)"""),
    CLI.Option("--phase-scaling-strength", "", """Sets the phase scaling strength. Corresponds to power
        values for positive, negative and triangular phase scaling
        type. (default: 4)"""),
    CLI.Option("--echoes", "-e", "Load only the specified echoes from disk (default: :)"; nargs=:many),
    CLI.Option("--no-mmap", "-N", """Has no effect: the inputs are always read into memory.
        Kept so that existing command lines still work."""; nargs=:none),
    CLI.Option("--no-phase-rescale", "", """Deactivate automatic rescaling of phase images. By
        default the input phase is rescaled to the range [-π;π]."""; nargs=:none),
    CLI.Option("--fix-ge-phase", "", """GE systems write corrupted phase output (slice jumps).
        This option fixes the phase problems."""; nargs=:none),
    CLI.Option("--writesteps", "", """Set to the path of a folder, if intermediate steps should
        be saved."""),
    CLI.Option("--verbose", "-v", "verbose output messages"; nargs=:none),
]

# The command line as given; `no_mmap` is also set when an input is compressed.
Base.@kwdef mutable struct ClearswiOptions
    magnitude::String = ""
    phase::String = ""
    output::String = "clearswi.nii"
    echo_times::Vector{String} = String[]
    mip_slices::String = "7"
    qsm::Bool = false
    qsm_input::String = ""
    qsm_mask::String = ""
    mag_combine::Vector{String} = ["SNR"]
    mag_sensitivity_correction::String = "on"
    mag_softplus_scaling::String = "on"
    unwrapping_algorithm::String = "laplacian"
    filter_size::Vector{String} = ["[4,4,0]"]
    phase_scaling_type::String = "tanh"
    phase_scaling_strength::String = "4"
    echoes::Vector{String} = [":"]
    no_mmap::Bool = false
    no_phase_rescale::Bool = false
    fix_ge_phase::Bool = false
    writesteps::String = ""
    verbose::Bool = false
end

_one(v::CLI.Values, key, default::String) = haskey(v, key) ? v[key][1] : default
_many(v::CLI.Values, key, default::Vector{String}) = haskey(v, key) ? v[key] : default
_flag(v::CLI.Values, key) = haskey(v, key)

function ClearswiOptions(v::CLI.Values)
    return ClearswiOptions(;
        magnitude = _one(v, "magnitude", ""),
        phase = _one(v, "phase", ""),
        output = _one(v, "output", "clearswi.nii"),
        echo_times = _many(v, "echo-times", String[]),
        mip_slices = _one(v, "mip-slices", "7"),
        qsm = _flag(v, "qsm"),
        qsm_input = _one(v, "qsm-input", ""),
        qsm_mask = _one(v, "qsm-mask", ""),
        mag_combine = _many(v, "mag-combine", ["SNR"]),
        mag_sensitivity_correction = _one(v, "mag-sensitivity-correction", "on"),
        mag_softplus_scaling = _one(v, "mag-softplus-scaling", "on"),
        unwrapping_algorithm = _one(v, "unwrapping-algorithm", "laplacian"),
        filter_size = _many(v, "filter-size", ["[4,4,0]"]),
        phase_scaling_type = _one(v, "phase-scaling-type", "tanh"),
        phase_scaling_strength = _one(v, "phase-scaling-strength", "4"),
        echoes = _many(v, "echoes", [":"]),
        no_mmap = _flag(v, "no-mmap"),
        no_phase_rescale = _flag(v, "no-phase-rescale"),
        fix_ge_phase = _flag(v, "fix-ge-phase"),
        writesteps = _one(v, "writesteps", ""),
        verbose = _flag(v, "verbose"),
    )
end

function getargs(args::AbstractVector, version)
    isempty(args) && (args = ["--help"])
    values = CLI.parse(CLI.Spec("clearswi", string(version), OPTIONS), args)
    values === nothing && return nothing
    return ClearswiOptions(values)
end

"""
    clearswi_main(args; version)

Runs clearswi with the command line `args` and returns the exit code. See
`clearswi_main(["--help"])` for the options.
"""
function clearswi_main(args; version=package_version(@__MODULE__))
    opts = try
        getargs(args, version)
    catch e
        e isa ArgumentError || rethrow()
        msg = e.msg
        msg isa String || (msg = "wrong arguments")
        CLI.print_stderr("clearswi: " * msg * "\nwrong argument formatting! See clearswi --help\n")
        return 1
    end
    opts === nothing && return 0 # --help or --version
    run_clearswi(opts, args, version)
    return 0
end

_nothing_if_empty(s::String) = isempty(s) ? nothing : s

function run_clearswi(opts::ClearswiOptions, args, version)
    writedir, filename = if endswith(opts.output, ".nii") || endswith(opts.output, ".nii.gz")
        dirname(opts.output), basename(opts.output)
    else
        opts.output, "clearswi"
    end
    if !isempty(opts.phase) && (endswith(opts.phase, ".gz") || endswith(opts.magnitude, ".gz"))
        opts.no_mmap = true
    end
    writesteps = _nothing_if_empty(opts.writesteps)

    mkpath(writedir)
    writesteps === nothing || mkpath(writesteps)

    isempty(opts.magnitude) && error("no magnitude image given (-m)")
    mag, hdr = loadmag(opts.magnitude)
    qsm_input = !isempty(opts.qsm_input)
    phase, phase_hdr = if qsm_input
        loadmag(opts.qsm_input) # use qsm instead of phase
    else
        isempty(opts.phase) && error("no phase image given (-p)")
        loadphase(opts.phase; rescale=!opts.no_phase_rescale, fix_ge=opts.fix_ge_phase)
    end
    phase_dims = Int(phase_hdr.dim[1])
    if opts.fix_ge_phase && writesteps !== nothing
        _in_dims(phase, phase_dims) do p
            savenii(p, "corrected_GE_phase", writesteps, hdr)
        end
    end
    neco = size(mag, 4)

    ## Echoes for unwrapping
    echoes = try
        _parse_echoes(opts.echoes, neco)
    catch y
        if isa(y, BoundsError)
            error("echoes=$(join(opts.echoes, " ")): specified echo out of range! Number of echoes is $neco")
        else
            error("echoes=$(join(opts.echoes, " ")) wrongly formatted!")
        end
    end
    opts.verbose && CLI.info("Echoes are " * _format_scalar_or_vector(echoes))

    TEs = _getTEs(opts, neco, echoes)
    opts.verbose && CLI.info("TEs are " * CLI.format(TEs))

    ## Error messages
    if 1 < length(echoes) && length(echoes) != length(TEs)
        error("Number of chosen echoes is $(length(echoes)) ($neco in .nii data), but $(length(TEs)) TEs were specified!")
    end

    filter_size = _parse_numbers(opts.filter_size)

    # Written here rather than before loading, so the record holds the resolved
    # echo times and echo count instead of only the raw arguments, and so the
    # citation list can depend on them.
    saveconfiguration(writedir, opts, filter_size, TEs, echoes, neco, args, version)

    phase_hp_sigma = _float_vector(filter_size)
    mip_slices = parse(Int, opts.mip_slices)
    phase_scaling_type = Symbol(opts.phase_scaling_type)
    phase_unwrap = Symbol(opts.unwrapping_algorithm)
    qsm = qsm_input ? :input : opts.qsm
    mag_softplus = if opts.mag_softplus_scaling == "on"
        true
    elseif opts.mag_softplus_scaling == "off"
        false
    else
        error("The setting for mag-softplus-scaling is not valid: $(opts.mag_softplus_scaling)")
    end
    phase_options = PhaseSettings(phase_unwrap, phase_hp_sigma, phase_scaling_type, writesteps, qsm)

    # A single echo is processed as 3D data; so is a QSM given as input once
    # echoes are selected. Otherwise the data keep the dimensions of the files.
    # Every branch below hands on arrays of static type.
    if length(echoes) == 1
        e = echoes[1]
        m = mag[:,:,:,e,1]
        p = qsm_input ? phase[:,:,:,1,1] : phase[:,:,:,e,1]
        _clearswi(opts, m, p, hdr, TEs, phase_options, mag_softplus, writedir, filename, mip_slices)
    elseif echoes == 1:neco
        m = _as4d(mag)
        if phase_dims <= 3
            _clearswi(opts, m, phase[:,:,:,1,1], hdr, TEs, phase_options, mag_softplus, writedir, filename, mip_slices)
        else
            _clearswi(opts, m, _as4d(phase), hdr, TEs, phase_options, mag_softplus, writedir, filename, mip_slices)
        end
    else
        opts.verbose && CLI.info("Selecting echoes " * CLI.format(echoes))
        m = mag[:,:,:,echoes,1]
        if qsm_input
            _clearswi(opts, m, phase[:,:,:,1,1], hdr, TEs, phase_options, mag_softplus, writedir, filename, mip_slices)
        else
            _clearswi(opts, m, phase[:,:,:,echoes,1], hdr, TEs, phase_options, mag_softplus, writedir, filename, mip_slices)
        end
    end
    return nothing
end

# The phase options that have one type whatever the command line
struct PhaseSettings
    phase_unwrap::Symbol
    phase_hp_sigma::Vector{Float64}
    phase_scaling_type::Symbol
    writesteps::Union{String,Nothing}
    qsm::Union{Bool,Symbol}
end

function _clearswi(opts, mag, phase, hdr, TEs, phase_options, mag_softplus, writedir, filename, mip_slices)
    data = Data(NoCopy(), mag, phase, hdr, TEs)
    swimag = _swimag(opts, data, mag_softplus, phase_options.writesteps)
    swiphase = _swiphase(opts, data, phase_options)
    _write_swi(swimag, swiphase, hdr, writedir, filename, mip_slices)
end

# Each part is Float32 or Float64, depending on the options. A function
# barrier rather than a broadcast in place, so that every combination is
# compiled with static types.
function _write_swi(swimag, swiphase, hdr, writedir, filename, mip_slices)
    swi = swimag .* swiphase
    mip = createIntensityProjection(swi, minimum, mip_slices)

    savenii(swi, filename, writedir, hdr)
    savenii(mip, "mip", writedir, hdr)
end

# The magnitude options, split into their types one at a time
function _swimag(opts, data, mag_softplus, writesteps)
    c = opts.mag_combine
    if c[1] == "SNR"
        return _swimag(opts, data, :SNR, mag_softplus, writesteps)
    elseif c[1] == "average"
        return _swimag(opts, data, :average, mag_softplus, writesteps)
    elseif c[1] == "echo"
        return _swimag(opts, data, :echo => parse(Int, last(c)), mag_softplus, writesteps)
    elseif c[1] == "SE"
        return _swimag(opts, data, :SE => parse(Float32, last(c)), mag_softplus, writesteps)
    end
    error("The setting for mag-combine is not valid: $(join(c, " "))")
end

function _swimag(opts, data, mag_combine, mag_softplus, writesteps)
    s = opts.mag_sensitivity_correction
    if s == "on"
        return getswimag(data, mag_combine, nothing, mag_softplus, writesteps)
    elseif s == "off"
        return getswimag(data, mag_combine, [1], mag_softplus, writesteps)
    elseif isfile(s)
        sensitivity = first(loadmag(s))
        return getswimag(data, mag_combine, sensitivity[:,:,:,1,1], mag_softplus, writesteps)
    end
    error("The setting for mag-sensitivity-correction is not valid: $s")
end

# The phase options, split into their types one at a time
function _swiphase(opts, data, o)
    strength = tryparse(Int, opts.phase_scaling_strength)
    if strength === nothing
        return _swiphase(opts, data, o, parse(Float32, opts.phase_scaling_strength))
    end
    return _swiphase(opts, data, o, strength)
end

function _swiphase(opts, data, o, strength)
    if !isempty(opts.qsm_mask)
        qsm_mask = first(loadmag(opts.qsm_mask))[:,:,:,1,1] .!= 0
        return getswiphase(data, PhaseOptions(o.phase_unwrap, o.phase_hp_sigma, o.phase_scaling_type, strength, o.writesteps, o.qsm, qsm_mask, nothing))
    end
    return getswiphase(data, PhaseOptions(o.phase_unwrap, o.phase_hp_sigma, o.phase_scaling_type, strength, o.writesteps, o.qsm, nothing, nothing))
end

_as4d(a::AbstractArray{<:Any,5}) = reshape(a, Base.front(size(a)))

# The data in the dimensions of the file
function _in_dims(f::F, a::Array{T,5}, nd) where {F,T}
    nd <= 3 && return f(a[:,:,:,1,1])
    nd == 4 && return f(_as4d(a))
    return f(a)
end

function _parse_echoes(strs, neco)
    parsed = MriResearchTools.parse_array(strs)
    echoes = if parsed isa Colon
        collect(1:neco)
    elseif parsed isa Int
        [parsed]
    elseif parsed isa Vector{Int}
        parsed
    else
        throw(ArgumentError("echoes must be integers"))
    end
    return (1:neco)[echoes] # throws BoundsError for an echo out of range
end

function _parse_numbers(strs)
    parsed = MriResearchTools.parse_array(strs)
    parsed isa Int && return [parsed]
    parsed isa Float64 && return [parsed]
    parsed isa Vector{Int} && return parsed
    parsed isa Vector{Float64} && return parsed
    throw(ArgumentError("\"$(join(strs, " "))\" is not a number or an array of numbers"))
end

_float_vector(v::AbstractVector{<:Real}) = Float64.(v)

function _getTEs(opts::ClearswiOptions, neco, echoes)
    if isempty(opts.echo_times)
        if neco == 1 || length(echoes) == 1
            return [1]
        else
            error("multi-echo data is used, but no echo times are given. Please specify the echo times using the -t option.")
        end
    end
    TEs = if opts.echo_times[1] == "epi"
        ones(neco) .* (length(opts.echo_times) > 1 ? parse(Float64, opts.echo_times[2]) : 1.0)
    else
        _parse_numbers(opts.echo_times)
    end
    if length(TEs) == neco
        TEs = TEs[echoes]
    end
    return TEs
end

# a single echo or echo time is recorded as the number itself
_format_scalar_or_vector(v::AbstractVector) = length(v) == 1 ? CLI.format(v[1]) : CLI.format(v)

function saveconfiguration(writedir, opts::ClearswiOptions, filter_size, TEs, echoes, neco, args, version)
    f = CLI.format
    settings = Dict{String,String}(
        "magnitude" => f(opts.magnitude),
        "phase" => f(opts.phase),
        "output" => f(opts.output),
        "echo-times" => f(opts.echo_times),
        "mip-slices" => f(opts.mip_slices),
        "qsm" => f(opts.qsm),
        "qsm-input" => f(opts.qsm_input),
        "qsm-mask" => f(opts.qsm_mask),
        "mag-combine" => f(opts.mag_combine),
        "mag-sensitivity-correction" => f(opts.mag_sensitivity_correction),
        "mag-softplus-scaling" => f(opts.mag_softplus_scaling),
        "unwrapping-algorithm" => f(opts.unwrapping_algorithm),
        "filter-size" => _format_scalar_or_vector(filter_size),
        "phase-scaling-type" => f(opts.phase_scaling_type),
        "phase-scaling-strength" => f(opts.phase_scaling_strength),
        "echoes" => f(opts.echoes),
        "no-mmap" => f(opts.no_mmap),
        "no-phase-rescale" => f(opts.no_phase_rescale),
        "fix-ge-phase" => f(opts.fix_ge_phase),
        "writesteps" => f(opts.writesteps),
        "verbose" => f(opts.verbose),
        "resolved-echo-times" => f(TEs),
        "resolved-echoes" => _format_scalar_or_vector(echoes),
        "number-of-echoes" => f(neco),
    )

    # Cite only what this run used.
    cite = [:clearswi]
    if opts.mag_sensitivity_correction == "on"
        push!(cite, :homogeneity)
    end
    if opts.unwrapping_algorithm == "romeo"
        push!(cite, :romeo)
    elseif contains(opts.unwrapping_algorithm, "laplacian")
        push!(cite, :laplacian)
    end
    packages = [CLEARSWI, MriResearchTools, MriResearchTools.ROMEO]
    if opts.qsm && isempty(opts.qsm_input)
        # The QSM path goes through QuantitativeSusceptibilityMappingTGV, which is
        # what the compiled app depends on. With --qsm-input the user supplies a
        # finished map, so nothing is cited for it here.
        append!(cite, [:tgv, :tgv_original])
        # qsm_romeo_B0 unwraps with ROMEO, and removes phase offsets with
        # MCPC-3D-S first when there is more than one echo.
        push!(cite, :romeo)
        if neco > 1
            push!(cite, :mcpc3ds)
        end
        # The QSM implementation is a weak dependency of MriResearchTools and
        # cannot be named here, so the module it registered its version for is
        # looked up.
        for m in keys(MriResearchTools.PACKAGE_VERSIONS)
            nameof(m) === :QuantitativeSusceptibilityMappingTGV && push!(packages, m)
        end
    end

    write_provenance(abspath(writedir), "clearswi";
        version, args, settings, cite,
        optional = [:julia],
        inputs = Pair{String,Union{Nothing,String}}["magnitude" => _nothing_if_empty(opts.magnitude), "phase" => _nothing_if_empty(opts.phase)],
        packages,
        describe = describe_input,
    )
end
