#!/usr/local/bin/julia
#title           :TverskyMatrix.jl
#description     :Calculate Tversky index similarity matrix from feature matrix
#author          :Carlos Vigil Vásquez
#date            :20220829
#version         :20220829a
#notes           :Requires ArgParse.jl & ProgressMeter.jl. `julia -t NUM_THREADS` for multithreading
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (carlos.vigil.v@gmail.com).
#license         :Permission to copy and modify is granted under the MIT license

using Base
using ArgParse
using DelimitedFiles
using ProgressMeter
using NamedArrays
using .Threads
include("named_array_helper.jl")

Base.setdiff(A::AbstractVector{Bool}, B::AbstractVector{Bool}) = @. A & !B

"""
    Tversky(X::AbstractVector{Bool}, Y::AbstractVector{Bool}, α::Number, β::Number)

Calculate the Tversky index between a pair of bit vectors.

# Arguments
- `M::AbstractMatrix{Bool}`: Bit matrix.
- `N::AbstractMatrix{Bool}`: Bit matrix.
- `α::AbstractFloat`: α coefficient.
- `β::AbstractFloat`: β coefficient.

# Examples
```jldoctest
julia> m = Vector{Bool}([1, 0, 0, 1, 0])
5-element Vector{Bool}:
 1
 0
 0
 1
 0

julia> n = Vector{Bool}([1, 0, 1, 1, 0])
5-element Vector{Bool}:
 1
 0
 1
 1
 0

julia> Tversky(m, n, 0.5, 0.5)    # Equivalent to Sørensen–Dice coefficient
0.8

julia> Tversky(m, n, 1.0, 1.0)    # Equivalent to Tanimoto coefficient
0.6666666666666666

julia> Tversky(m, n, 1.0, 0.0)    # "Superstructure-likeness" measure
1.0

julia> Tversky(m, n, 0.0, 1.0)    # "Substructure-likeness" measure
0.6666666666666666
```

# Extended help
Tversky index is equivalent to other similarity coefficients for some special cases:
- For `α = β = 1.0`, Tversky index is equivalent to Tanimoto coefficient
- For `α = β = 0.5`, Tversky index is equivalent to Sørensen–Dice coefficient
- For `α = 1.0` & `β = 0.0`, Tversky index corresponds to a "superstucture-likeness" measure
- For `α = 0.0` & `β = 1.0`, Tversky index corresponds to a "substucture-likeness" measure

Tversky index does not meet the criteria for a similarity metric.
"""
function Tversky(X::T, Y::T, α::U, β::U) where {T<:AbstractVector{Bool},U<:Number}
    # Check if coefficients are greater or equal to zero
    @assert α ≥ 0 "α must be positive!"
    @assert β ≥ 0 "β must be positive!"

    # Check if both matrices have the same number of columns
    @assert length(X) == length(Y)

    return sum(X .& Y) / (sum(X .& Y) + α * sum(setdiff(X, Y)) + β * sum(setdiff(Y, X)))
end

"""
    Tversky(M::AbstractMatrix{Bool}, N::AbstractMatrix{Bool}, α::Number, β::Number)

Calculate the Tversky index between a pair of bit matrices.

# Arguments
- `M::AbstractMatrix{Bool}`: Bit matrix.
- `N::AbstractMatrix{Bool}`: Bit matrix.
- `α::AbstractFloat`: α coefficient.
- `β::AbstractFloat`: β coefficient.

# Examples
```jldoctest
julia> M = Matrix{Bool}([0 1 0 0 1; 1 0 1 0 0; 1 0 0 0 1; 0 0 1 0 0; 0 0 0 1 1])
5×5 Matrix{Bool}:
 0  1  0  0  1
 1  0  1  0  0
 1  0  0  0  1
 0  0  1  0  0
 0  0  0  1  1

julia> Tversky(M, M, 0.5, 0.5)
5×5 Matrix{Float64}:
 1.0  0.0       0.5  0.0       0.5
 0.0  1.0       0.5  0.666667  0.0
 0.5  0.5       1.0  0.0       0.5
 0.0  0.666667  0.0  1.0       0.0
 0.5  0.0       0.5  0.0       1.0

julia> Tversky(M, M, 1.0, 0.0)
5×5 Matrix{Float64}:
 1.0  0.0  0.5  0.0  0.5
 0.0  1.0  0.5  0.5  0.0
 0.5  0.5  1.0  0.0  0.5
 0.0  1.0  0.0  1.0  0.0
 0.5  0.0  0.5  0.0  1.0
```

# Extended help
Tversky index is equivalent to other similarity coefficients for some special cases:
- For `α = β = 1.0`, Tversky index is equivalent to Tanimoto coefficient
- For `α = β = 0.5`, Tversky index is equivalent to Sørensen–Dice coefficient
- For `α = 1.0` & `β = 0.0`, Tversky index corresponds to a "superstucture-likeness" measure
- For `α = 0.0` & `β = 1.0`, Tversky index corresponds to a "substucture-likeness" measure

Tversky index does not meet the criteria for a similarity metric.
"""
function Tversky(M::T, N::T, α::U, β::U) where {T<:AbstractMatrix{Bool},U<:Number}
    # Check if coefficients are greater or equal to 0
    @assert α ≥ 0 "α must be positive!"
    @assert β ≥ 0 "β must be positive!"

    # Check if both matrices have the same number of columns
    Mₘ, Mₙ = size(M)
    Nₘ, Nₙ = size(N)
    @assert Mₙ == Nₙ "Matrices should have the same amount of columns!"

    # Create measurement matrix
    MN = zeros(Mₘ, Nₘ)
    pbar = Progress(Mₘ; showspeed=true)
    @threads for i in 1:Mₘ
        for j in 1:Nₘ
            @inbounds MN[i, j] = Tversky(M[i, :], N[j, :], α, β)
        end
        next!(pbar)
    end

    return MN
end

function main(args)
    # Argument parsing {{{
    configs = ArgParseSettings()

    add_arg_group!(configs, "I/O option:")
    @add_arg_table! configs begin
        "--MD", "-i"
        arg_type = String
        action = :store_arg
        help = "M × D matrix"
        required = true
        "--ND"
        arg_type = String
        action = :store_arg
        help = "N × D matrix"
        required = false
        "--MN", "-o"
        arg_type = String
        action = :store_arg
        help = "M × N similarity matrix"
        required = true
        "--named"
        action = :store_true
        help = "Matrix has index/names in first column"
        "--delimiter", "-d"
        arg_type = Char
        action = :store_arg
        help = "Delimiter character between bits"
        required = false
        default = ' '
    end

    add_arg_group!(configs, "Tversky index options (defaults to Tanimoto coefficient):")
    @add_arg_table! configs begin
        "--alpha"
        arg_type = Float64
        action = :store_arg
        help = "α coefficient"
        default = 1.0
        "--beta"
        arg_type = Float64
        action = :store_arg
        help = "β coefficient"
        default = 1.0
    end

    args = parse_args(args, configs)
    args["ND"] = args["ND"] === nothing ? args["MD"] : args["ND"]
    # }}}

    # Load matrices and get number of drugs, compounds, substructures and targets
    println("Calculating Tversky index (α = $(args["alpha"]), β = $(args["beta"])) matrix from file $(args["ND"])")
    MD = parse_matrix(readdlm(args["MD"], args["delimiter"], String), true, false; type=Bool)
    ND = parse_matrix(readdlm(args["ND"], args["delimiter"], String), true, false; type=Bool)

    # Calculate similarity matrix and save
    MN = Tversky(MD, ND, args["alpha"], args["beta"])
    if args["named"]
        MN = NamedArray(MN)
        setnames!(MN, names(MD, 1), 1)
        setnames!(MN, names(ND, 1), 2)
    end
    writedlm(args["MN"], MN, args["delimiter"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
