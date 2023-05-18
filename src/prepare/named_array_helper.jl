using NamedArrays
import DelimitedFiles.writedlm

writedlm(io::IO, x::NamedMatrix{T} where {T}) = writedlm(io, vcat(["" names(x, 2)...], hcat(names(x, 1), x)))
writedlm(io::IO, x::NamedMatrix{T} where {T}, delimiter::Char) = writedlm(io, vcat(["" names(x, 2)...], hcat(names(x, 1), x)), delimiter)
writedlm(io::AbstractString, x::NamedMatrix{T} where {T}, delimiter::Char) = writedlm(io, vcat(["" names(x, 2)...], hcat(names(x, 1), x)), delimiter)

"""
    parse_matrix(M::AbstractMatrix, rows::Bool, cols::Bool; type::Type)

Convert matrix with row/column names to NamedMatrix (assumes names are in firt row/column).

# Arguments
- `M::AbstractMatrix` : Matrix to parse
- `rows::Bool` : Matrix has row names (default = false)
- `cols::Bool` : Matrix has column names (default = false)
- `type::Type` : Type of matrix values (default = Any)
"""
function parse_matrix(M::AbstractMatrix, rows::Bool=false, cols::Bool=false; type::Type=Any)
    # Extract values from matrix
    c_idx = rows ? 2 : 1
    r_idx = cols ? 2 : 1
    values = parse.(type, M[r_idx:end, c_idx:end])

    # Extract dimensions names
    row_names = rows ? [i for i in String.(M[r_idx:end, 1])] : ["R#$i" for i in 1:size(values, 1)]
    col_names = cols ? [i for i in String.(M[1, c_idx:end])] : ["C#$i" for i in 1:size(values, 2)]

    namedM = NamedArray(values, (row_names, col_names))
    namedM = namedM[sort(row_names), sort(col_names)]

    return namedM
end

function main()
    A = [1 1 0 1; 0 1 1 0; 1 0 1 1]
    B = ["A" 1 1 0 1 1; "B" 1 0 0 1 0; "C" 1 0 0 1 1]
    C = ["F1" "F2" "F3" "F4"; 1 1 0 1; 0 1 1 0; 1 0 1 1]
    D = ["" "F1" "F2" "F3" "F4"; "A" 1 1 0 1; "B" 0 1 1 0; "C" 1 0 1 1]

    parse_matrix(A, false, false)
    parse_matrix(B, true, false)
    parse_matrix(C, false, true)
    parse_matrix(D, true, true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
