
__precompile__()

module CopyNum

using Helpers, Printf, Statistics, DelimitedFiles
export CopyNumber, Region, Segment, writetsv
export mean_decimate, median_decimate, segment
export cn_to_logratio, logratio_to_cn, read_igv, write_igv, read_segments

mutable struct CopyNumber
	logratio::Matrix{Float32}
	logratio_ci::Matrix{Float32}   # 95% confidence interval
	median_af::Matrix{Float32}
	call::Matrix{Int16}
	hetz_snps::Matrix{Int32}    # Number of hetz SNPs used for median AF calc
	num_probes::Vector{Int32}
	sample::Vector{String}
	sample_noise::Vector{Float32}
	gene::Vector{String}
	chromosome::Vector{String}
	position::Matrix{Int32}
end

struct Region
	chromosome::String
	start::Int32
	stop::Int32
end

mutable struct Segment
	chromosome::String
	start::Int32
	stop::Int32
	logratio::Float32
	baf::Float32
	logratio_var::Float32
	baf_var::Float32
end

cn_to_logratio(cn::Real, cancer_frac::Real, ploidy::Int) =
	log2(cancer_frac * (cn / ploidy - 1) + 1)
function logratio_to_cn(lr::Real, cancer_frac::Real, ploidy::Int)
	if !(0 <= cancer_frac <= 1)
		error("Cancer fraction must be between 0 and 1.")
	end
	return ploidy * (1 + (2^lr - 1) / cancer_frac)
end

Base.size(cn::CopyNumber, args...) = Base.size(cn.logratio, args...)
Base.getindex(cn::CopyNumber, rows, cols) = CopyNumber(
	cn.logratio[rows, cols], cn.logratio_ci[rows, cols],
	cn.median_af[rows, cols], cn.call[rows, cols],
	cn.hetz_snps[rows, cols], cn.num_probes[rows],
	cn.sample[cols], cn.sample_noise[cols], 
	cn.gene[rows], cn.chromosome[rows], cn.position[rows, :])

#function Base.getindex(cn::CopyNumber, gene_name::AbstractString, cols)
#	row = get(cn.gene_name_to_row, gene_name, 0)
#	if row == 0; error("Gene '$(gene_name)' not found."); end
#	return cn[row, cols]
#end

function CopyNumber(genes::AbstractVector, samples::AbstractVector)
	S = length(samples); G = length(genes)
	#gene_name_to_row = Dict(g => row for (row, g) in enumerate(genes))
	return CopyNumber(zeros(Float32, G, S), zeros(Float32, G, S),
		fill(Float32(NaN), G, S),
		zeros(Int16, G, S), zeros(Int32, G, S), zeros(Int32, G),
		samples, fill(Float32(NaN), S), genes, fill("", G), zeros(Int32, G, 2))
end

function CopyNumber(tsv_path::AbstractString)
	d = readdlm(tsv_path, '\t'); headers = d[1, :][:];
	genes = d[2:end, 1][:]; samples = headers[5:end];
	cn = CopyNumber(genes, samples)
	for g in 1:length(genes)
		cn.chromosome[g] = d[1+g, 2]
		cn.position[g, 1] = d[1+g, 3]
		cn.position[g, 2] = d[1+g, 4]
		for s in 1:length(samples)
			fields = split(d[1+g, 4+s], ':')
			cn.call[g, s] = parse(Int, fields[1])
			cn.logratio[g, s] = parse(Float32, fields[2])
			cn.median_af[g, s] = parse(Float32, fields[3])
			cn.hetz_snps[g, s] = parse(Int, fields[4])
		end
	end
	return cn
end

function writetsv(out::IO, cn::CopyNumber)
	println("Sample\tGene\tChrom\tStart\tEnd\tCoverage logratio\tCoverage logratio 95% CI\tMedian heterozygous SNP allele fraction\tNumber of heterozygous SNPs")
	for s in 1:length(cn.sample), g in 1:length(cn.gene)
		@printf("%s\t%s\t%s\t%d\t%d", cn.sample[s], cn.gene[g],
			cn.chromosome[g], cn.position[g, 1], cn.position[g, 2])
		@printf("\t%.3f\t%.3f\t%.3f\t%d\n", cn.logratio[g, s],
			cn.logratio_ci[g, s], cn.median_af[g, s], cn.hetz_snps[g, s])
	end
end
writetsv(out_path::AbstractString, cn::CopyNumber; kwargs...) =
	open(fd -> writetsv(fd, cn; kwargs...), expanduser(out_path), "w")

function mean_decimate(values::Vector, fold::Integer)
	if fold <= 1; return copy(values); end
	if isempty(values); return zeros(eltype(values), 0); end
	starts = 1:fold:length(values)
	decimated = zeros(Float64, length(starts))
	for k in 1:length(starts)-1     # Last window handled as special case
		decimated[k] = mean(values[(0:fold-1) .+ starts[k]])
	end
	decimated[length(starts)] = mean(values[starts[end]:end])
	return decimated
end

function mean_decimate(chromosome::Vector{String}, position::Vector, value::Vector, fold::Int)
	decimated_chr = fill("", 0)
	decimated_pos = zeros(Int, 0)
	decimated_value = zeros(eltype(value), 0)

	for run in runs(chromosome)
		pos = round.(Int, mean_decimate(Float64.(position[run]), fold))
		append!(decimated_chr, fill(chromosome[run[1]], length(pos)))
		append!(decimated_pos, pos)
		append!(decimated_value, mean_decimate(value[run], fold))
	end
	return decimated_chr, decimated_pos, decimated_value
end

function median_decimate(values::Vector, fold::Integer)
	if fold <= 1; return copy(values); end
	if isempty(values); return zeros(eltype(values), 0); end
	starts = 1:fold:length(values)
	decimated = zeros(eltype(values), length(starts))
	for k in 1:length(starts)-1     # Last window handled as special case
		decimated[k] = median(values[(0:fold-1) .+ starts[k]])
	end
	decimated[length(starts)] = median(values[starts[end]:end])
	return decimated
end

function median_decimate(chromosome::Vector{String}, position::Vector, value::Vector, fold::Int)
	decimated_chr = fill("", 0)
	decimated_pos = zeros(Int, 0)
	decimated_value = zeros(eltype(value), 0)

	for run in runs(chromosome)
		pos = round.(Int, median_decimate(Float64.(position[run]), fold))
		append!(decimated_chr, fill(chromosome[run[1]], length(pos)))
		append!(decimated_pos, pos)
		append!(decimated_value, median_decimate(value[run], fold))
	end
	return decimated_chr, decimated_pos, decimated_value
end

# Old circular binary segmentation Z(). The new ones is non-circular.
#Z(S, i, j) = ((S[j] - S[i]) / (j - i) - (S[end] - S[j] + S[i]) / (length(S) - j + i)) / sqrt(1 / (j - i) + 1 / (length(S) - j + i))

Z(S, i, j) = (S[i] / i - (S[j] - S[i]) / (j - i)) / sqrt(1 / i + 1 / (j - i))

# Generated using: round.(10 .* 1.5.^(0:17))
const window_sizes = [10, 15, 22, 34, 51, 76, 114, 171, 256, 384, 577, 865, 1297, 1946, 2919, 4379, 6568, 9853, 300_000_000];

function find_breakpoint(data::Vector)
	N = length(data)
	S = cumsum(data)
	best_z = 0
	best_split = 1
	for split in 1:N - 1
		for wsize in window_sizes
			i = split; j = min(split + wsize, N)
			z = abs(Z(S, i, j))            # Max			
			if z > best_z
				best_z = z
				best_split = split
			end
			if j == length(data); break; end
		end
	end
	return (best_split, best_z)
end

function segment(chr, position, logratios; threshold=1)
	chr_ranges = Dict(c => findfirst(chr .== c):findlast(chr .== c) for c in unique(chr))

	S = size(logratios, 2)
	noise_levels = [median(abs.(diff(logratios[:, s]))) for s in 1:S]

	segments = Region[]
	for (chr, range) in chr_ranges
		# Recursively split segments until we are done
		segments_todo = [range]
		while !isempty(segments_todo)
			segment = pop!(segments_todo)
			best_split = 1; best_z = 0
			for s in 1:S
				split, z = find_breakpoint(logratios[segment, s])
				z /= (40 * noise_levels[s])
				if z > best_z
					best_split = split
					best_z = z
				end
			end

			if best_z >= threshold
				push!(segments_todo, segment[1:best_split])
				push!(segments_todo, segment[best_split+1:end])
			else
				push!(segments, Region(chr, segment[1], segment[end]))
			end
		end
	end
	sort!(segments, by=(s -> s.start))

	# Calculate median logratios and variances for each segment
	regions = segments
	segments = Matrix{Segment}(undef, length(regions), S)
	for (r, reg) in enumerate(regions), s in 1:S
		valid_lr = filter(x -> !isnan(x), logratios[reg.start:reg.stop, s])
		segments[r, s] = Segment(reg.chromosome, position[reg.start],
			position[reg.stop], isempty(valid_lr) ? NaN : median(valid_lr),
			NaN, isempty(valid_lr) ? NaN : var(valid_lr), NaN)
	end
	return segments
end

function read_igv(igv_path::AbstractString)
	d = readtsv(igv_path)
	chr = String.(d[2:end, 1])
	pos = [round(Int, (d[r, 2] + d[r, 3]) / 2) for r in 2:size(d, 1)]
	value = Float32.(d[2:end, 5])
	return (chr, pos, value)
end

function write_igv(path::String, samples::Vector, chromosome::Vector, position::Vector, values::Array)
	out = open(expanduser(path), "w")
	write(out, "CHROM\tSTART\tEND\tFEATURE\t$(join(samples, '\t'))\n")
	for k in 1:length(chromosome)
		@printf(out, "%s\t%d\t%d\t", chromosome[k], position[k], position[k]+1)
		for v in values[k, :]; @printf(out, "\t%.3f", v); end
		write(out, '\n')
	end
	close(out)
end

function read_segments(seg_file; min_size=1e6, skip_xy=false)
	d = readtsv(seg_file)
	segments = Segment[]
	for k in 2:size(d, 1)
		chr = String(d[k, 1])
		start = Float64(d[k, 2])
		stop = Float64(d[k, 3])
		if skip_xy && (chr == "chrX" || chr == "chrY"); continue; end
		if stop - start < min_size; continue; end
		
		if ':' in d[k, 4]
			parts = split(d[k, 4], ':')
			logratio = parse(Float32, parts[1])
			saf = parse(Float32, parts[2])
			logratio_var = parse(Float32, parts[3])
			saf_var = parse(Float32, parts[4])
			push!(segments, Segment(chr, start, stop, logratio, saf, logratio_var, saf_var))
		else
			logratio = d[k, 4]
			push!(segments, Segment(chr, start, stop, logratio, NaN, NaN, NaN))
		end
	end
	return segments
end

end