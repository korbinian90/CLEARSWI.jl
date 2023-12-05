function getargs(args::AbstractVector, version)
    if isempty(args)
        args = ["--help"]
    end
    s = ArgParseSettings(
        exc_handler=exception_handler,
        add_version=true,
        version,
        )
    @add_arg_table! s begin
        "--magnitude", "-m"
            help = "The magnitude image (single or multi-echo)"
        "--phase", "-p"
            help = "The phase image (single or multi-echo)"
        "--output", "-o"
            help = "The output path or filename"
            default = "clearswi.nii"
        "--echo-times", "-t"
            help = """The echo times are required for multi-echo datasets 
                specified in array or range syntax (eg. "[1.5,3.0]" or 
                "3.5:3.5:14")."""
            nargs = '+'
        "--qsm"
            help = """When activated uses RTS QSM for phase weighting.
            """
            action = :store_true
        "--mag-combine"
            help = """SNR | average | echo <n> | SE <te>.
                Magnitude combination algorithm. echo <n> selects a specific
                echo; SE <te> simulates a single echo scan of the given echo
                time."""
            default = ["SNR"]
            nargs = '+'
        "--mag-sensitivity-correction"
            help = """ <filename> | on | off.
                Use the CLEAR-SWI sensitivity correction. Alternatively, a
                sensitivity map can be read from a file"""
            default = "on"
        "--mag-softplus-scaling"
            help = """on | off.
                Set softplus scaling of the magnitude"""
            default = "on"
        "--unwrapping-algorithm"
            help = """laplacian | romeo | laplacianslice"""
            default = "laplacian"
        "--filter-size"
            help = """Size for the high-pass phase filter in voxels. Can be
                given as <x> <y> <z> or in array syntax (e.g. [2.2,3.1,0],
                which is effectively a 2D filter)."""
            nargs = '+'
            default = ["[4,4,0]"]
        "--phase-scaling-type"
            help = """tanh | negativetanh | positive | negative | triangular
                Select the type of phase scaling. positive or negative with a
                strength of 3-6 is used in standard SWI."""
            default = "tanh"
        "--phase-scaling-strength"
            help = """Sets the phase scaling strength. Corresponds to power
                values for positive, negative and triangular phase scaling
                type."""
            default = "4"
        "--echoes", "-e"
            help = "Load only the specified echoes from disk"
            default = [":"]
            nargs = '+'
        "--no-mmap", "-N"
            help = """Deactivate memory mapping. Memory mapping might cause
                problems on network storage"""
            action = :store_true
        "--no-phase-rescale"
            help = """Deactivate automatic rescaling of phase images. By
                default the input phase is rescaled to the range [-π;π]."""
            action = :store_true
        "--writesteps"
            help = """Set to the path of a folder, if intermediate steps should
                be saved."""
            default = nothing
        "--verbose", "-v"
            help = "verbose output messages"
            action = :store_true
    end
    return parse_args(args, s)
end

function exception_handler(settings::ArgParseSettings, err, err_code::Int=1)
    if err == ArgParseError("too many arguments")
        println(stderr,
            """wrong argument formatting!"""
        )
    end
    ArgParse.default_handler(settings, err, err_code)
end

function getechoes(settings, neco)
    echoes = eval(Meta.parse(join(settings["echoes"], " ")))
    if echoes isa Int
        echoes = [echoes]
    elseif echoes isa Matrix
        echoes = echoes[:]
    end
    echoes = (1:neco)[echoes] # expands ":"
    if (length(echoes) == 1) echoes = echoes[1] end
    return echoes
end

function getTEs(settings, neco, echoes)
    if isempty(settings["echo-times"])
        if neco == 1 || length(echoes) == 1
            return [1]
        else
            error("multi-echo data is used, but no echo times are given. Please specify the echo times using the -t option.")
        end
    end
    TEs = if settings["echo-times"][1] == "epi"
        ones(neco) .* if length(settings["echo-times"]) > 1; parse(Float64, settings["echo-times"][2]) else 1 end
    else
        eval(Meta.parse(join(settings["echo-times"], " ")))
    end
    if TEs isa Matrix
        TEs = TEs[:]
    end
    if length(TEs) == neco
        TEs = TEs[echoes]
    end
    if !(TEs isa AbstractVector)
        TEs = [TEs]
    end
    return TEs
end

function saveconfiguration(writedir, settings, args, version)
    writedir = abspath(writedir)
    open(joinpath(writedir, "settings_clearswi.txt"), "w") do io
        for (fname, val) in settings
            if !(typeof(val) <: AbstractArray)
                println(io, "$fname: " * string(val))
            end
        end
        println(io, """Arguments: $(join(args, " "))""")
        println(io, "ClearswiApp version: $version")
    end
end
