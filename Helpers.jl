
__precompile__()

module Helpers

using DelimitedFiles

export readtsv, findone, hierarchical_order, runs, info, warn
export readtsv_dict

info(msg::String) = println(stderr, msg)
warn(msg::String) = println(stderr, "WARNING: $msg")

# Allow logical operators on boolean arrays without the dot operator
Base.:!(a::AbstractArray{Bool}) = .!a
Base.:&(a::AbstractArray{Bool}, b::AbstractArray{Bool}) = a .& b
Base.:|(a::AbstractArray{Bool}, b::AbstractArray{Bool}) = a .| b

Base.findfirst(s::AbstractString, c::Char) = Base.findfirst(isequal(c), s)

function findone(predicate::Function, collection::AbstractVector)
	idx = 0
	for k in 1:length(collection)
		if predicate(collection[k])
			if idx > 0
				error("At least two elements $(collection[idx]) and $(collection[k]) fulfill the predicate.")
			end
			idx = k
		end
	end
	return idx == 0 ? nothing : idx
end

function findone(collection::AbstractVector, elem)
	idx = 0
	for k in 1:length(collection)
		if collection[k] == elem
			if idx > 0
				error("Element '$elem' is found in at least two positions $idx and $k.")
			end
			idx = k
		end
	end
	return idx == 0 ? nothing : idx
end

function Base.only(predicate::Function, collection::AbstractVector)
	idx = findone(predicate, collection)
	if idx == 0; error("No elements fulfill the predicate."); end
	return collection[idx]
end

# Allow writing replace(some_strings, "pattern", "replacement")
Base.replace(strings::AbstractArray, pattern, sub) =
	map(s -> replace(s, pattern => sub), strings)
Base.replace(s::AbstractString, pattern, sub) = replace(s, pattern => sub)

# These allow in() and the in-operator to be used to test if a string
# contains a pattern (where pattern can be a string or regular expression).
Base.in(pattern::String, s::AbstractString) = occursin(pattern, s)
Base.in(pattern::SubString{String}, s::AbstractString) = occursin(pattern, s)
Base.in(r::Regex, s::AbstractString) = occursin(r, s)

readtsv(tsv_file::IO; text=false) = readdlm(tsv_file, '\t', text ? String : Any)
readtsv(cmd::Base.AbstractCmd; kwargs...) =
	open(f -> readtsv(f; kwargs...), cmd)
readtsv(tsv_path::AbstractString; kwargs...) =
	tsv_path == "-" ? readtsv(STDIN; kwargs...) : open(f -> readtsv(f; kwargs...), expanduser(tsv_path))

function hierarchical_order(args...)
	for a in 2:length(args); @assert(length(args[a]) == length(args[1])); end
	function lt(a, b)
		for c in 1:length(args)
			if args[c][a] < args[c][b]; return true; end
			if args[c][a] > args[c][b]; return false; end
		end
		return false
	end
	return sortperm(1:length(args[1]), lt=lt)
end

function runs(values::AbstractVector)
	runs = Vector{UnitRange{Int}}()
	start = 1
	for k in 2:length(values)
		if values[k] != values[k-1]
			push!(runs, start:(k-1))
			start = k
		end
	end
	push!(runs, start:length(values))
	return runs
end

end