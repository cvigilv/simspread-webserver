#!/usr/local/bin/julia

using ArgParse
using CUDA
using DelimitedFiles
using LinearAlgebra
using NamedArrays
using ProgressMeter
using SimSpread

function main(args)
    # Argument parsing {{{
    configs = ArgParseSettings()

    add_arg_group!(configs, "I/O:")
    @add_arg_table! configs begin
        "--dt-train"
        help = "Training drug-target adjacency matrix"
        required = true
        action = :store_arg
        arg_type = String
        "--dd-train"
        help = "Training drug-drug similarity matrix"
        required = true
        action = :store_arg
        arg_type = String
        "--dd-query"
        help = "Testing drug-drug similarity matrix"
        required = true
        action = :store_arg
        arg_type = String
        "--output-file", "-o"
        help = "Predicted test drug-target interactions"
        required = true
        action = :store_arg
        arg_type = String
    end

    add_arg_group!(configs, "wl-SimSpread parameters:")
    @add_arg_table! configs begin
        "--weighted", "-w"
        help = "Similarity matrix featurization weighting scheme"
        action = :store_true
        "--cutoff", "-c"
        help = "Similarity cutoff (α)"
        action = :store_arg
        arg_type = Float64
        default = 0.5
    end

    add_arg_group!(configs, "Miscellaneous:")
    @add_arg_table! configs begin
        "--gpu"
        help = "GPU acceleration"
        action = :store_true
        "--gpu-id"
        help = "GPU to use"
        action = :store_arg
        arg_type = Int64
        default = 0
    end

    parsed_args = parse_args(args, configs)

    # Store arguments to variables
    weighted = parsed_args["weighted"]
    output = parsed_args["output-file"]
    α = parsed_args["cutoff"]
    # }}}

    # Load matrices to memory
    DD₀ = read_namedmatrix(parsed_args["dd-train"])
    DD₁ = read_namedmatrix(parsed_args["dd-query"])
    DT₀ = read_namedmatrix(parsed_args["dt-train"])
    DT₁ = NamedArray(zeros(1, size(DT₀, 2)), [(names(DD₁, 1)), (names(DT₀, 2))])

    # Predict interactions
    DF₀ = featurize(DD₀, α, weighted)
    DF₁ = featurize(DD₁, α, weighted)
    G = construct((DT₀, DT₁), (DF₀, DF₁))
    R = predict(G, DT₁; GPU=CUDA.functional())
    open(output, "w+") do f
        for (st, score) in enamerate(R)
            s,t = st
            write(f, join([s,t,score], ' ') * '\n')
        end
    end

end

main(ARGS)
