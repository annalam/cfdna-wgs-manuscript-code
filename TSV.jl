
__precompile__()

module TSV

using DelimitedFiles, Dates

export read_tsv, readtsv_dict

# Helpers for reading tab-delimited files
mutable struct DelimitedFile
	headers::Vector{String}
	header_index::Dict{String, Int}
	data::Matrix{Any}
end
Base.size(f::DelimitedFile, args...) = Base.size(f.data, args...)
Base.lastindex(f::DelimitedFile, dim::Integer) = Base.size(f.data, dim)
Base.getindex(f::DelimitedFile, rows, col::AbstractString) = 
	Base.getindex(f.data, rows, get(f.header_index, col, 0))
Base.getindex(f::DelimitedFile, rows, col_regex::Regex) =
	Base.getindex(f.data, rows, col_regex .∈ f.headers)
Base.getindex(f::DelimitedFile, rows, cols) = Base.getindex(f.data, rows, cols)
Base.show(io::IO, f::DelimitedFile) = Base.show(io, f.data)
#TODO? Base.setindex(...)

function read_tsv(tsv_file::IO; header=true, text=false)
	d = readdlm(tsv_file, '\t', text ? String : Any)
	headers = header ? map(x -> "$x", d[1, :][:]) : []
	header_index = Dict(zip(headers, 1:length(headers)))
	return DelimitedFile(headers, header_index, header ? d[2:end, :] : d)
end
read_tsv(cmd::Base.AbstractCmd; kwargs...) =
	open(f -> read_tsv(f; kwargs...), cmd)
read_tsv(tsv_path::AbstractString; kwargs...) =
	tsv_path == "-" ? read_tsv(STDIN; kwargs...) : open(f -> read_tsv(f; kwargs...), expanduser(tsv_path))
function readtsv_dict(tsv_path::AbstractString; dim="cols")
	d=readdlm(tsv_path, '\t')
	if dim .== "rows"
		keynames = d[:,1]
		d = permutedims(d[:,2:end], (2,1))
	elseif dim .== "cols"
		keynames = string.(d[1,:])
		d = d[2:end,:]
	end
	for key in unique(keynames)
		dups = findall(keynames.==key)
		if length(dups) > 1
			for (i, di) in enumerate(dups)
				keynames[di] = keynames[di]*"("*string(i)*")"
			end
		end
	end
	Dict(keynames[ki]=>d[:,ki] for ki in 1:length(keynames))
end

end