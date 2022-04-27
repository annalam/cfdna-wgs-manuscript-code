
__precompile__()

module CopyNum

using Helpers, Printf, Statistics, DelimitedFiles, Serialization, Tracks
using Distributions: Normal, quantile, Binomial, Poisson

export CopyNumber, Segment, writetsv
export mean_decimate, median_decimate, segment
export cn_to_logratio, logratio_to_cn
export simulate_hsaf_deviation, hsaf_deviation_from_table,  deviated_hsaf_table, expected_hsaf
export true_hsaf_from_table
export expected_maf, mutant_allele_cn
export nearest_cn_state, cn_state, cn_states, mut_segment_cn_states
export read_igv, write_igv, read_segments, read_genome_track

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

mutable struct Segment
	chromosome::String
	start::Int32
	stop::Int32
	logratio::Float32
	hsaf::Float32
	logratio_var::Float32
	hsaf_stdev::Float32
	depth::Float32		# Median sequencing depth at HSAF sites
end

cn_to_logratio(cn::Real, cancer_frac::Real, ploidy::Integer) =
	log2(cancer_frac * (cn / ploidy - 1) + 1)

# This function assumes that coverage logratios have been normalized such
# that logratio = 0 corresponds to CN = normal ploidy.
function logratio_to_cn(lr::Real, cancer_frac::Real, ploidy::Integer)
	if !(0 <= cancer_frac <= 1)
		error("Cancer fraction must be between 0 and 1.")
	end
	return ploidy * (1 + (2^lr - 1) / cancer_frac)
end

# This function calculates the expected heterozygous SNP allele fraction for
# a given allele-specific copy number state.
global deviated_hsaf_table = nothing
function hsaf_deviation_from_table(true_hsaf, depth)
	global deviated_hsaf_table
    if deviated_hsaf_table == nothing
    	data_dir = "$(@__DIR__)/../data"
    	deviated_hsaf_table = deserialize("$data_dir/deviated_hsaf_table.jls")["deviated_hsaf_table"]
    end
	depth_i = max(1, round.(Int, depth))
	depth_i = min(size(deviated_hsaf_table, 1), depth_i)
    hsaf_i = round.(Int, abs(true_hsaf-0.5)*1000)+1
    hsaf = deviated_hsaf_table[depth_i, hsaf_i]::Float64
end
function simulate_hsaf_deviation(true_hsaf, depth; iterations=10_000)
    if typeof(true_hsaf) != Float64
        true_hsaf = Float64(true_hsaf)
    end
    deviation = zeros(iterations)
    d = 0
    step = ceil(Int, iterations / 100)
    for k in 1:step:iterations
    	# We do Float64(depth) here because recent versions of the
    	# Distributions.jl package do not allow a Float32 depth.
        dpt = rand(Poisson(Float64(depth)))
        dist = Binomial(dpt, true_hsaf)
        for _ in 1:min(step, iterations - k + 1)
        	d += 1
            deviation[d] = 0.5 + abs(0.5 - rand(dist) / dpt)
        end
    end
    hsaf = median(deviation)
    if true_hsaf < 0.5
        hsaf = 1 - hsaf
    end
    return hsaf
end
function expected_hsaf(cancer_frac::Real, b_allele_cn::Real, total_cn::Real; depth=0, iterations=10_000, use_table=false)
	@assert(b_allele_cn <= total_cn)
	true_hsaf = (cancer_frac * (b_allele_cn - 1) + 1) / (cancer_frac * (total_cn - 2) + 2)
	if depth == 0
		return true_hsaf
	elseif use_table
		return hsaf_deviation_from_table(true_hsaf, depth)
	else
		# Run a simulation to estimate the expected HSAF 
		return simulate_hsaf_deviation(true_hsaf, depth; iterations=10_000)
	end
end

# This function returns estimate for the true HSAF of a given observed HSAF and depth
global true_hsaf_table = nothing
function true_hsaf_from_table(obs_hsaf, depth)
    if true_hsaf_table == nothing
    	data_dir = "$(@__DIR__)/../data"
    	true_hsaf_table = deserialize("$data_dir/true_hsaf_table.jls")["true_hsaf_table"]
    end
	depth_i = max(1, round.(Int, depth))
	depth_i = min(size(true_hsaf_table, 1), depth_i)
    hsaf_i = round.(Int, abs(obs_hsaf-0.5)*1000)+1
    hsaf = true_hsaf_table[depth_i, hsaf_i]::Float64
end

# This function calculates the expected somatic mutation allele fraction for
# a truncal somatic mutation with a given mutant allele copy number in a genomic
# region with a given total copy number.
expected_maf(cancer_frac::Real, alt_cn::Real, tot_cn::Real, normal_cn::Real) =
	(cancer_frac*alt_cn) / (cancer_frac*tot_cn + (1-cancer_frac)*normal_cn)

mutant_allele_cn(maf::Real, cancer_frac::Real, cn::Real, normal_cn::Int) =
	maf * (cn + normal_cn * (1 / cancer_frac - 1))

# Returns a (total copy number, major allele copy number) tuple.
# Assumes that logratios have been normalized such that logratio = 0
# corresponds to CN = 2. Only applicable to autosomes.
function nearest_cn_state(logratio::Real, hsaf::Real, cancer_frac::Real)
	cn = round(Int, logratio_to_cn(logratio, cancer_frac, 2))
	if cn <= 0; return (0, 0); end
	possible_major_cn = ceil(Int, cn / 2):cn
	best = Int(0); best_delta = Inf
	for major_cn in possible_major_cn
		exp_hsaf = expected_hsaf(cancer_frac, major_cn, cn)
		delta = abs(hsaf - exp_hsaf)
		if delta < best_delta
			best = major_cn
			best_delta = delta
		end
	end
	return (cn, best)
end

# This function checks if there is a truncal allele-specific copy number state
# that matches well with the coverage logratio and heterozygous SNP allele
# fractions of a given genomic segment.
#
# Returns a (total copy number, major allele copy number) tuple.
# Assumes that logratios have been normalized such that logratio = 0
# corresponds to CN = 2. Only applicable to autosomes.
function cn_state(cancer_frac::Real, median_depth::Real, lr::Real, 
	lr_var::Real, hsaf::Real, hsaf_std::Real; 
	q_threshold=NaN, lr_q_threshold=0.15, hsaf_q_threshold=0.15, 
	separation_threshold=NaN, cn_separation_threshold=0.25, hsaf_separation_threshold=0.25, 
	hsaf_iterations=10_000, use_hsaf_table=false)
	
	if !isnan(q_threshold)
		lr_q_threshold=q_threshold
		hsaf_q_threshold=q_threshold
	end
	if !isnan(separation_threshold)
		cn_separation_threshold=separation_threshold
		hsaf_separation_threshold=separation_threshold
	end
	T = Union{NTuple{4, Nothing}, Tuple{Int, Int, Float64, Float64}}

    # Find the total copy number state that best matches the coverage logratio
    raw_cn = logratio_to_cn(lr, cancer_frac, 2)
    cn::Int = round(Int, raw_cn)
    if cn < 0; cn = 0; end # This line is causing some weird type problems. The above ::Int helped a little

    # Check that the measured logratio is close enough to the mathematically
	# expected logratio, allowing for some error in the measurement
	expected_lr = cn_to_logratio(cn, cancer_frac, 2)
	lr_dist = Normal(lr, sqrt(lr_var))
    if !(quantile(lr_dist, lr_q_threshold) < expected_lr < quantile(lr_dist, 1 - lr_q_threshold))
        return T((nothing, nothing, nothing, nothing))
    end

	# Check that the next closest CN level is not too close
	if !(cn - cn_separation_threshold <= raw_cn <= cn + cn_separation_threshold)
        return T((nothing, nothing, nothing, nothing))
	end

	# nearest allele-specific level
	possible_major_cn = Int.(ceil(cn / 2):cn)
	possible_hsaf = [expected_hsaf(cancer_frac, major_cn, cn, depth=median_depth, 
						iterations=hsaf_iterations, use_table=use_hsaf_table)
						for major_cn in possible_major_cn]
    best = argmin(abs.(possible_hsaf .- hsaf))
    exp_hsaf = possible_hsaf[best]
	major_cn = possible_major_cn[best]

	# Check that the next closest CN level is not too close
	if hsaf < exp_hsaf && best!=1
		if hsaf < exp_hsaf - hsaf_separation_threshold * (exp_hsaf - possible_hsaf[best-1])
			return T((nothing, nothing, nothing, nothing))
		end
	elseif hsaf > exp_hsaf && best != length(possible_hsaf)
		if hsaf > exp_hsaf + hsaf_separation_threshold * (possible_hsaf[best+1] - exp_hsaf)
			return T((nothing, nothing, nothing, nothing))
		end
	end

    # check that the level is inside hsaf quantile thresholds
    hsaf_dist = Normal(hsaf, hsaf_std)
    if !(quantile(hsaf_dist, hsaf_q_threshold) < exp_hsaf < quantile(hsaf_dist, 1 - hsaf_q_threshold))
        return T((nothing, nothing, nothing, nothing))
    end

	return T((cn, major_cn, abs(expected_lr-lr), abs(exp_hsaf-hsaf)))
end

# Returns allele specific copy number statuses for all segments in a vector
function cn_states(cancer_frac::Real, diploid_level::Real, segments::AbstractVector{Segment};
			max_cn=Inf, q_threshold=0.15, separation_threshold=0.25,
            use_hsaf_table=false, hsaf_iterations=10_000)
	tots = Vector{Union{Nothing, Int64}}(undef, length(segments))
	majs = Vector{Union{Nothing, Int64}}(undef, length(segments))
    for (s, seg) in enumerate(segments)
        if isnan(seg.logratio_var) || isnan(seg.hsaf_stdev); continue; end
        tot, maj = cn_state(cancer_frac, seg.depth, seg.logratio-diploid_level, 
			seg.logratio_var, seg.hsaf, seg.hsaf_stdev, q_threshold=q_threshold, 
			separation_threshold=separation_threshold, hsaf_iterations=hsaf_iterations, 
			use_hsaf_table=use_hsaf_table)
        if tot==nothing || tot > max_cn; continue; end
		tots[s] = tot
		majs[s] = maj
    end
	return (tot_cn=tots, maj_cn=majs)
end

function mut_segment_cn_states(mut_chr::Vector, mut_pos::Vector, segments::Vector{Segment}, 
								seg_tot_cns::Vector, seg_maj_cns::Vector)
	m_tot_cns = Vector{Union{Int64, Nothing}}(undef, length(mut_chr));
	m_maj_cns = Vector{Union{Int64, Nothing}}(undef, length(mut_chr));
	seg = 1
	for mut in 1:length(mut_chr), seg in 1:length(segments)
		if mut_chr[mut] != segments[seg].chromosome; continue; end
		if !(segments[seg].start < mut_pos[mut] < segments[seg].stop); continue; end
		m_tot_cns[mut] = seg_tot_cns[seg]
		m_maj_cns[mut] = seg_maj_cns[seg]
	end
	return (tot_cn=m_tot_cns, maj_cn=m_maj_cns)
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

	if eltype(values) <: Integer
		# For integer inputs, if the window size is even, we may need to
		# round the output values to integers as well.
		for k in 1:length(starts)-1     # Last window handled as special case
			decimated[k] = round(eltype(values),
				median(values[(0:fold-1) .+ starts[k]]))
		end
		decimated[length(starts)] = round(eltype(values),
			median(values[starts[end]:end]))
	else
		# For floating point inputs, we don't need to round the outputs
		for k in 1:length(starts)-1     # Last window handled as special case
			decimated[k] = median(values[(0:fold-1) .+ starts[k]])
		end
		decimated[length(starts)] = median(values[starts[end]:end])
	end

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

median_decimate(gt::GenomeTrack, fold::Int) =
	GenomeTrack(median_decimate(gt.chromosome, gt.position, gt.value, fold)...)

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
			NaN, isempty(valid_lr) ? NaN : var(valid_lr), NaN, NaN)
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
		
		logratio = Float64(d[k, 4])
		hsaf = Float64(d[k, 6])
		logratio_var = Float64(d[k, 5])
		hsaf_stdev = Float64(d[k, 7])
		depth = size(d, 2) >= 8 ? Float64(d[k, 8]) : NaN
		push!(segments, Segment(chr, start, stop, logratio, hsaf, logratio_var, hsaf_stdev, depth))
	end
	return segments
end

end