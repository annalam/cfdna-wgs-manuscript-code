
__precompile__()

module Tracks

using Helpers

export Region, GenomeTrack, read_genome_track

struct Region
	chromosome::String
	start::Int32
	stop::Int32
end

struct GenomeTrack{T<:Real}
	chromosome::Vector{String}
	position::Vector{Int}
	value::Vector{T}
	function GenomeTrack{T}(chr::Vector{String}, pos::Vector{Int}, val::Vector{T}) where T <: Real
		@assert length(chr) == length(pos) == length(val)
		new(chr, pos, val)
	end
end

GenomeTrack(chr, pos, val::AbstractVector{T}) where T <: Real =
	GenomeTrack{T}(string.(chr), Int.(pos), collect(val))

Base.length(gt::GenomeTrack) = length(gt.position)
#Base.size(gpv::GenomeTrack) = size(gpv.value)
#Base.size(gpv::GenomeTrack, d::Int) = size(gpv.value, d)
#Base.ndims(gpv::GenomeTrack) = ndims(gpv.value)
Base.lastindex(gt::GenomeTrack) = length(gt.position)

Base.getindex(gt::GenomeTrack, i) = Base.getindex(gt.value, i)


# Create a GenomeTrack from a input with no headers and three columns:
# chr, pos and value.
# If a file includes headers and/or columns are in different order, user can do for example:
# 	read_genome_track(`cut -f 1,2,4 input.tsv`)
# 	read_genome_track(pipeline(`tail +2 input.tsv`, `cut -f 1,2,4`))
function read_genome_track(tsv::Union{IO,Base.AbstractCmd,AbstractString}; header=false)
	tsv = readtsv(tsv)
	first_row = header ? 2 : 1
	GenomeTrack(tsv[first_row:end, 1], tsv[first_row:end, 2], Float64.(tsv[first_row:end, 3]))
end

# Create a GenomeTrack from multiple inputs akin to the previous method's input
# Inputs must have same number of lines with identical chr and pos columns (identicality not checked!)
# Values are stored in a matrix with columns corresponding to the different inputs
# function read_genome_track(tsvs::AbstractVector{T} where T <: Union{IO,Base.AbstractCmd,AbstractString}; header=false)
# 	tsv1 = readtsv(tsvs[1])
# 	if header; tsv1 = tsv1[2:end,:]; end
# 	chr, pos, lr1 = collect(eachcol(tsv1))
# 	lr = zeros(typeof(lr1[1]), length(lr1), length(tsvs))
# 	lr[:,1] .= lr1
# 	for i in 2:length(tsvs)
# 		lr_s = readtsv(tsvs[i])[:,end]
# 		if header; lr_s = lr_s[2:end]; end
# 		lr[:,i] .= lr_s
# 	end
# 	return GenomeTrack(chr, pos, lr)
# end

# function read_regions(bed_path)
# 	chromosomes = Dict{Any, Regions}()
# 	d = readdlm(expanduser(bed_path), '\t')
# 	for chr in unique(d[:, 1])
# 		pos = int(d[d[:, 1] .== chr, 2:3]); pos[:, 1] += 1
# 		strand = size(d, 2) >= 6 ? d[:, 6] .== "+" : trues(pos)
# 		chromosomes[replace("$(chr)", "chr", "")] = Regions(pos, strand)
# 	end
# 	return chromosomes
# end

# near(pos::Integer, regions::Regions, radius::Integer) = 
# 	any(k -> regions.position[k, 1] - radius <= pos <= regions.position[k, 2] + radius, 1:length(regions))

# near(chromosomes::Array, positions::Array, regions::Dict{Any, Regions},
# 	radius::Integer) = map(k -> near(positions[k], regions[chromosomes[k]], radius), 1:length(chromosomes))

# near_upstream(pos::Integer, regions::Regions, radius::Integer) = 
# 	any(k -> abs(regions.position[k, regions.strand[k] ? 1 : 2] - pos) <= radius, 1:length(regions))

# near_upstream(chromosomes::Array, positions::Array, regions::Dict{Any, Regions}, radius::Integer) = map(k -> near_upstream(positions[k], regions[chromosomes[k]], radius), 1:length(chromosomes))

# max_near(chromosomes::Array, positions::Array, tracks::Dict{Any, TrackFixed}, radius::Integer) = convert(Array{Float32}, [maximum(tracks[chromosomes[k]][(-radius:radius) + positions[k]]) for k in 1:length(chromosomes)])



end