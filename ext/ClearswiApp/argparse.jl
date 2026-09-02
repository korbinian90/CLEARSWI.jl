function getargs(args::AbstractVector, version)
    if isempty(args)
        args = ["--help"]
    end
    s = ArgParseSettings(;
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
        "--mip-slices", "-s"
            help = "The number of slices in the MIP image"
            default = "7"
        "--qsm"
            help = """When activated uses TGV QSM for phase weighting.
            """
            action = :store_true
        "--qsm-input"
            help = """Give pre-calculated QSM instead of phase as input. Not fine-tuned!
                Please adjust phase filter to something like:
                "--filter-size [20,20,0] --phase-phase_scaling_strength 0.04" """
            default = nothing
        "--qsm-mask"
            help = """The mask used for QSM. Use a custom mask, if the qsm_mask.nii 
                is not good for your data."""
            default = nothing
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
        "--fix-ge-phase"
            help = """GE systems write corrupted phase output (slice jumps).
                This option fixes the phase problems."""
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
    echoes = MriResearchTools.ROMEO.parse_array(settings["echoes"])
    if echoes isa Int
        echoes = [echoes]
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
        MriResearchTools.ROMEO.parse_array(settings["echo-times"])
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
    # Cite only what this run used.
    cite = [:clearswi]
    if settings["mag-sensitivity-correction"] == "on"
        push!(cite, :homogeneity)
    end
    if settings["unwrapping-algorithm"] == "romeo"
        push!(cite, :romeo)
    elseif contains(settings["unwrapping-algorithm"], "laplacian")
        push!(cite, :laplacian)
    end
    if settings["qsm"] === true && isnothing(settings["qsm-input"])
        # The QSM path goes through QuantitativeSusceptibilityMappingTGV, which is
        # what the compiled app depends on. It used to cite the RTS dipole
        # inversion of QSM.jl, which this tool does not run. With --qsm-input the
        # user supplies a finished map, so nothing is cited for it here.
        append!(cite, [:tgv, :tgv_original])
        # qsm_romeo_B0 unwraps with ROMEO, and removes phase offsets with
        # MCPC-3D-S first when there is more than one echo.
        push!(cite, :romeo)
        if get(settings, "number-of-echoes", 1) > 1
            push!(cite, :mcpc3ds)
        end
    end

    # The QSM implementation is a weak dependency of MriResearchTools and cannot
    # be named here, so look the loaded module up when it ran.
    packages = Any[CLEARSWI, MriResearchTools, MriResearchTools.ROMEO]
    if :tgv in cite
        tgv = _loaded_module("QuantitativeSusceptibilityMappingTGV",
                             "bd393529-335a-4aed-902f-5de61cc7ff49")
        tgv === nothing || push!(packages, tgv)
    end

    write_provenance(writedir, "clearswi";
        version, args, settings, cite,
        optional = [:julia],
        inputs = ["magnitude" => settings["magnitude"], "phase" => settings["phase"]],
        packages,
        describe = describe_input,
    )
end

_loaded_module(name, uuid) =
    get(Base.loaded_modules, Base.PkgId(Base.UUID(uuid), name), nothing)
