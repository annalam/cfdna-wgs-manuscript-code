
__precompile__()

module Subclones

export find_truncal_model, truncal_model_plot, bayesian_clustering
export assign_mutation_clusters, clust_macn_to_cf, simulate_ccf
export mutation_genome_positions_plot
export solve_phylogeny, phylotree_plot
export cn_deconvolution, cn_deconvolution_plot
export clonal_expansion_plot
export CN_PALETTE

using Printf, Random, StatsBase
using Helpers, Plot2, DelimitedFiles, Statistics, Tracks, CopyNum
using Suppressor
using Distributions: Binomial, Normal, Categorical, pdf
using Interpolations: LinearInterpolation, SteffenMonotonicInterpolation, interpolate
using Combinatorics: combinations, permutations, partitions

const CN_PALETTE = Dict(nothing=>RGB(0, 255, 0), "subclonal"=>RGB(0, 255, 0), "bad"=>RGB(255, 0, 255),
						  0=>RGB(0, 0, 255), 1=>RGB(120, 120, 255), 2=>RGB(100), 3=>RGB(255, 180, 180),
						  4=>RGB(255, 120, 120), 5=>RGB(255, 0, 0), 6=>RGB(180, 0, 0), 7=>RGB(120, 0, 0), 
						  8=>RGB(60, 0, 0), 9=>RGB(40, 0, 0), 10=>RGB(20, 0, 0), 11=>RGB(10, 0, 0), 
						  12=>RGB(0), 13=>RGB(0), 14=>RGB(0)
)

# Human hg38 reference genome
const CHR_NAMES = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"];
const CHR_SIZES = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569];
const CHR_STARTS = cumsum(CHR_SIZES) .- CHR_SIZES;

function binomial_kde(grid, success_fractions::Vector, trials::Vector)
	density = zeros(length(grid))
	for i in 1:length(success_fractions)
		f = success_fractions[i]
		t = trials[i]
		if f==0 || t==0; continue; end
		p = pdf.(Binomial(t, f), 0:t)
		interp = LinearInterpolation((0:t)./t, p)
		density += interp.(grid)
	end
	density = density ./ length(success_fractions) .* 100
	return density
end


# TRUNCAL MODELS
logrange(start, stop, len) = 2 .^(range(log2(start), log2(stop), length=len))
function bincounts(x, bin_ceils)
	bcounts = zeros(eltype(bin_ceils), length(bin_ceils))
	x = sort!(x)
	b = 1
	for xx in x[searchsortedfirst(x, 2*bin_ceils[1]-bin_ceils[2]):end]
		while xx>bin_ceils[b]
			b+=1
		end
		bcounts[b] += 1
	end
	return bcounts
end
cn_logratio_to_cancerfrac(cn, lr) = 2*(2^lr-1)/(cn-2)
function convert_ploidy(from, to, cf, dl)
	dl_change = cn_to_logratio(2*from/to, cf, 2)
	cf_new = cn_logratio_to_cancerfrac(2*to/from, -dl_change)
	return cf_new, dl_change + dl
end
function find_truncal_model(segments::AbstractVector{T} where T <: Segment; n_ploidies=3,
							init_min_hsaf_prop=0.33, init_lr_bins=20, init_lr_peak_q=0.33, 
							frac_weight=1, dist_weight=1, impossible_hsaf_weight=50,
							q_threshold=NaN, lr_q_threshold=0.15, hsaf_q_threshold=0.15, 
							separation_threshold=NaN, cn_separation_threshold=0.25, hsaf_separation_threshold=0.25, 
							hsaf_dist_fac=1, max_0_copy_regions=50e6,
							require_minor_allele_0_segment=true)

	function model_error(cf::Float64, dl::Float64; prev_top_error=Inf,
						seg_len_norm=seg_len_norm, seg_len=seg_len, segments=segments)

		tots = zeros(Int, length(segments))
		majs = zeros(Int, length(segments))
		lr_dist = zeros(length(segments))
		hsaf_dist = zeros(length(segments))
		for (s, seg) in enumerate(segments)
			tot, maj, lrd, hsafd = cn_state(cf, seg.depth, seg.logratio.-dl, 
									seg.logratio_var, seg.hsaf, seg.hsaf_stdev,
									use_hsaf_table=true, q_threshold=q_threshold, 
									lr_q_threshold=lr_q_threshold, 
									hsaf_q_threshold=hsaf_q_threshold,
									separation_threshold=separation_threshold, 
									cn_separation_threshold=cn_separation_threshold, 
									hsaf_separation_threshold=hsaf_separation_threshold)
			tots[s], majs[s], lr_dist[s], hsaf_dist[s] = tot==nothing ? (-1, -1, NaN, NaN) : (tot, maj, lrd, hsafd)
		end

		if require_minor_allele_0_segment
			if !any(t!=-1 && t==m for (t,m) in zip(tots, majs)); return NaN; end # at least one segment with minor allele cn of 0
		end

		truncal = tots.!=-1
		if length(unique(zip(tots[truncal], majs[truncal]))) < 2; return NaN; end # at least two unique cn states
		assigned_frac_penalty = length(tots) / sum(truncal)
		error = frac_weight * assigned_frac_penalty
		# println(assigned_frac_penalty)
		if error >= prev_top_error; return NaN; end
	
		lr_dist = lr_dist[truncal]
		hsaf_dist = hsaf_dist[truncal]./cf
		assign_distance_penalty = sum(seg_len_norm[truncal] .* sqrt.(lr_dist .^2 .+ (hsaf_dist_fac .* hsaf_dist) .^2))
		error += dist_weight * assign_distance_penalty
		# println(assign_distance_penalty)
		if error >= prev_top_error; return NaN; end
	
		raw_cns = max.(0.0, logratio_to_cn.([s.logratio for s=segments].-dl, cf, 2))
		max_hsafs::Vector{Float64} = broadcast(raw_cns, [s.depth for s=segments]) do cn, depth
			expected_hsaf.(cf, cn, cn, depth=depth, use_table=true)
		end
		max_hsaf_ratios = max_hsafs ./ [s.hsaf for s=segments]
		imposs = max_hsaf_ratios .< 1
		impossible_hsaf_penalty = sum((1 .- max_hsaf_ratios[imposs]).^2 .* seg_len_norm[imposs])
		error += impossible_hsaf_weight * impossible_hsaf_penalty
		# println(impossible_hsaf_penalty)
		if error >= prev_top_error; return NaN; end
	
		return error
	end

	function find_model(cf_start::Float64, cf_radius::Float64, cf_step::Float64,
						dl_start::Float64, dl_radius::Float64, dl_step::Float64)
		if isnan(cf_step)
			cf_range = range(0.0, 0.5, length=60) .+ logrange(0.01, 0.5, 60)
		else
			cf_range = max(cf_start-cf_radius, 0.001):cf_step:min(cf_start+cf_radius, 1.0)
		end
		dl_range = dl_start-dl_radius:dl_step:dl_start+dl_radius
		top_cf, top_dl, top_err = NaN, NaN, Inf
		for cf in cf_range, dl in dl_range
			err = model_error(cf, dl, prev_top_error=top_err)
			if isnan(err); continue; end
			top_cf, top_dl, top_err = cf, dl, err
		end
		return top_cf, top_dl, top_err
	end

	function tune_model(cf::Float64, dl::Float64)
		cf_t, dl_t, err_t = cf, dl, NaN
		for x in 3 .^(2:-1:0)
			radius = round(0.001*2*x, digits=3)
			step = round(0.001*x, digits=3)
			cf_t, dl_t, err_t = find_model(cf_t, radius, step, dl_t, radius, step)
			if isnan(cf_t)
				return cf, dl, model_error(cf, dl)
			end
		end
		return cf_t, dl_t, err_t
	end

	function zero_copy_regions(cf::Float64, dl::Float64; segments=segments, 
							   raw_copy_thresh=0.9, seg_len=seg_len)
		min_one_copy_lr = cn_to_logratio(raw_copy_thresh, cf, 2) + dl
		zero_copy_regs = 0
		for (i, s) in enumerate(segments) 
			if s.logratio < min_one_copy_lr
				zero_copy_regs += seg_len[i]
			end
		end
		return zero_copy_regs
	end

	segments = [s for s=segments if !isnan(s.logratio_var) && 
					!isnan(s.hsaf_stdev) && !isnan(s.depth)]
	seg_len = [s.stop-s.start for s=segments]
	seg_len_norm = max.(5e6, seg_len) ./ max(5e6, maximum(seg_len))

	min_hsaf, max_hsaf = extrema([s.hsaf for s=segments])
	min_hsaf_thresh = min_hsaf + (max_hsaf-min_hsaf)*init_min_hsaf_prop
	min_hsaf_segs = [s.hsaf for s=segments] .<= min_hsaf_thresh
	diploid_candidates = min_hsaf_segs & (seg_len .> 2e6)
	min_candidate_lr, max_candidate_lr = Float64.(extrema([s.logratio for s=segments[diploid_candidates]]))
	bins = range(min_candidate_lr, max_candidate_lr, length=init_lr_bins)[2:end]
	bcounts = bincounts([s.logratio for s=segments[diploid_candidates]], bins)
	bincount_thresh = sort(unique(bcounts))[end-floor(Int, (length(unique(bcounts))-1)*init_lr_peak_q)]
	lr_peaks = round.(collect(bins)[bcounts.>=bincount_thresh] .- (bins[2]-bins[1])/2, digits=2)

	dl_radius = round((max_candidate_lr-min_candidate_lr) / init_lr_bins / 2, digits=5)
	dl_step = dl_radius / 3

	init_models = [find_model(0.0, 1.0, NaN, dl, dl_radius, dl_step) for dl=lr_peaks];
	# cf_start, cf_radius, cf_step, dl_start, dl_radius, dl_step = 0.0, 1.0, NaN, lr_peaks[1], dl_radius, dl_step

	cf1, dl1, err1 = init_models[argmin([m[3] for m=init_models])]
	if isnan(cf1)
		return [(cf1, dl1, err1)]
	end
	cf1, dl1, err1 = tune_model(cf1, dl1)
	
	models = NamedTuple{(:cancer_fraction, :diploid_level, :error), Tuple{Float64, Float64, Float64}}[]

	if zero_copy_regions(cf1, dl1) <= max_0_copy_regions
		push!(models, (cancer_fraction=cf1, diploid_level=dl1, error=err1))
	end
	ploidy_change_penalty = 0.1
	for (i, ploidy_change) in enumerate(4:2:2*n_ploidies)
		cf, dl = convert_ploidy(2, ploidy_change, cf1, dl1)
		cf, dl, err = tune_model(cf, dl)
		if zero_copy_regions(cf, dl) > max_0_copy_regions; continue; end
		err *= (1+ploidy_change_penalty)^i
		push!(models, (cancer_fraction=cf, diploid_level=dl, error=err))
	end
	models = models[sortperm([m[3] for m=models])]

	return models
end

function segment_cn_plot(cf, dl, segments::Vector{Segment})
	seg_tot_cns, _ = cn_states(cf, dl, segments, max_cn=12, use_hsaf_table=true)

	# draw segemets
	seg_tot = [logratio_to_cn(s.logratio - dl, cf, 2) for s=segments] # raw segment CN values
	seg_hsaf = [s.hsaf for s=segments]
	point_size = [1+((s.stop-s.start) / 1e7) for s=segments]
	color = [CN_PALETTE[cn] for cn in seg_tot_cns]
	scatter_plot(seg_tot, seg_hsaf, color=color, size=point_size)

	# draw expected positions of allele specific cn statuses
	tots = 0:min(ceil(Int, maximum(seg_tot)), 14)
	k = median(s.depth / 2^(s.logratio + 1 - dl) for s=segments if !isnan(s.depth))
	tot, hsaf = zeros(0), zeros(0)
	for t=tots, maj=ceil(Int, t/2):t
		push!(tot, t)
		push!(hsaf, expected_hsaf(cf, maj, t, depth=round(Int, k*(cf*t + (1-cf)*2))))
	end
	color = [lighter_color(CN_PALETTE[cn]) for cn in tot]
	scatter_plot(tot, hsaf, border_color=color, fill_color="none", size=150)

	ylabel("Hetz SNP allele fraction")
	xlabel("Copy number")
	grid(b=false)
	xlim(NaN, min(get_xlim(2), 12))
end
function maf_kdes(cancer_fraction, tot_r, alt_r, tot_cn; min_muts=0)
	figure() do
		order = sortperm(tot_cn)
		idx_maf = [(order[r], alt_r[order[r]] ./ tot_r[order[r]])
				   for r=runs(tot_cn[order]) if length(r)>=min_muts]
		xmax = length(idx_maf)==0 ? NaN : 0
		for (_, maf) in idx_maf
			mafmax = maximum(maf)
			if mafmax+0.1>xmax
				xmax=mafmax+0.1 
			end
		end
		for (i, maf) in idx_maf
			t_cn = tot_cn[i[1]]
			x = 0:0.001:1
			y = binomial_kde(x, maf, tot_r[i])
			y = y ./ maximum(y) .+ get_ylim()[2]
			line_plot(x, y, color=CN_PALETTE[t_cn], line_width=2)
			for a_cn in 1:t_cn
				l = expected_maf(cancer_fraction, a_cn, t_cn, 2)
				d = y[findfirst(x .> l)]
				label = "$(a_cn)/$t_cn"
				label_dir = -1 + 2 * (d-minimum(y) < 0.5)
				tick_end = d + label_dir * 0.2
				v_align = label_dir==1 ? "bottom" : "top"
				line_plot(fill(l ,2), [d, tick_end], color=CN_PALETTE[t_cn], line_width=2)
				if l < xmax
					text(l, tick_end + label_dir * 0.05, label, size=8, v_align=v_align)
				end
			end
			legend = "Copy number $t_cn\n$(length(i)) mutations"
			text(xmax, get_ylim()[2]-0.1, legend, size=8, 
					v_align="top", h_align="right", color=CN_PALETTE[t_cn])
		end
		grid(b=false)
		xlim(0, xmax)
		ylim(0.9, NaN)
		yticks([])
		ylabel("Density")
		xlabel("Somatic mutation allele fraction")
	end
end
function truncal_model_plot(cancer_fraction::Float64, diploid_level::Real,
							coverage_logratio::GenomeTrack, hetz_snp_fraction::GenomeTrack,
							m_tot::Vector, m_alt::Vector, m_tot_cns::Vector, segments::Vector{Segment},
							seg_tot_cns::Vector, seg_maj_cns::Vector;
							title_prefix=nothing, max_cn=9, min_visualized_muts=100,
							adaptive_axis_limits=true, hsaf_visualization_quantile=1.0, # the y-axis limits of the visulaization
							lr_visualization_quantile=1.0, # "middle" quantile: e.g. 0.9 -> 0.05 to 0.95; only used if adaptive_axis_limits=false
							show_quantile_thresholds=false # segment CN assignment thresholds as vertical lines
							)

	# set max cn to max number in color palette
	max_cn = min(max_cn, maximum(filter(k->typeof(k)==Int, keys(CN_PALETTE))))
	for i in 1:length(seg_tot_cns)
		if seg_tot_cns[i]==nothing; continue;end
		if seg_tot_cns[i]>max_cn
			seg_maj_cns[i] = nothing
			seg_tot_cns[i] = nothing
		end
	end

	# All observed total copy numbers
	obs_cn_statuses = vcat([[maj tot] for (maj, tot) in unique(zip(seg_maj_cns, seg_tot_cns))
							if tot!=nothing]...)

	
	figure(size=(9,8)) do

		# Coverage logratio scatter plot
		subplot(1, 4)
			if adaptive_axis_limits && size(obs_cn_statuses, 1)!=0
				min_cn = min(1, minimum(obs_cn_statuses[:,2]))
				min_cn_l = cn_to_logratio(min_cn, cancer_fraction, 2)
				min2_cn_l = cn_to_logratio(min_cn+1, cancer_fraction, 2)
				ymin = min_cn_l - (min2_cn_l - min_cn_l)/2 + diploid_level
				max_cn = max(4, maximum(obs_cn_statuses[:,2]))
				if max_cn!=0
					max_cn_l = cn_to_logratio(max_cn, cancer_fraction, 2)
					max2_cn_l = cn_to_logratio(max_cn-1, cancer_fraction, 2)
					ymax = max_cn_l + (max_cn_l - max2_cn_l)  + diploid_level
				end
			else
				ymin = quantile(coverage_logratio.value, 0.5 - lr_visualization_quantile/2)
				ymax = quantile(coverage_logratio.value, 0.5 + lr_visualization_quantile/2)
			end

			# Plot horizontal lines for expected levels
			if adaptive_axis_limits
				tots = 1:4
			else
				totmin = round(Int, logratio_to_cn(minimum(s.logratio for s=segments), cancer_fraction, 2))
				totmax = round(Int, logratio_to_cn(maximum(s.logratio for s=segments), cancer_fraction, 2))
				tots = max(0, totmin):min(max_cn, totmax)
			end
			if size(obs_cn_statuses, 1)!=0
				tots = unique(vcat(tots, obs_cn_statuses[:,2]))
			end
			for tot_cn in tots
				l = cn_to_logratio(tot_cn, cancer_fraction, 2) + diploid_level
				if l < ymin; continue; end
				line_plot([0, CHR_STARTS[23]], fill(l, 2), color=CN_PALETTE[tot_cn], line_width=0.5)
				text(CHR_STARTS[25]+3e7, l, tot_cn, h_align="left", size=7, color=CN_PALETTE[tot_cn])
			end

			# Plot segment medians and thresholds
			for i in 1:length(segments)
				s = segments[i]
				if CHR_STARTS[findone(CHR_NAMES, s.chromosome)] in CHR_STARTS[23:25]; continue; end;
				start, stop = CHR_STARTS[findone(CHR_NAMES, s.chromosome)] .+ [s.start, s.stop]
				lr_mean, var = s.logratio, s.logratio_var
				color=CN_PALETTE[seg_tot_cns[i]]
				line_plot([start, stop], fill(lr_mean, 2), color=color, line_width=2)
				if show_quantile_thresholds && isfinite(var)
					@assert(!isnan(lr_q_thresh) && !isnan(hsaf_q_thresh),
							 "You must give thresholds to show them")
					@assert lr_q_thresh<0.5 "quantile thresholds must be between 0 and 5"
					# Horizontal line showing segment cn assignment thresholds
					lq, uq = quantile(Normal(lr_mean, sqrt(var)), [lr_q_thresh, 1-lr_q_thresh])
					line_plot(fill((start+stop)/2, 2), [lq, uq], color=color, line_width=1.2)
				end
			end

			# Scatter plot of coverage logratio values
			genome_scatter_plot(coverage_logratio.chromosome,
								coverage_logratio.position,
								coverage_logratio.value)

			ylim(ymin, ymax)
			ylabel("Coverage logratio")

		# HSAF scatter plot
		subplot(2, 4)
			ymax = round(quantile(Float64.(hetz_snp_fraction.value[:, 1]), hsaf_visualization_quantile) + 0.05, digits=1)

			# Plot horizontal lines for expected levels
			k = median(s.depth / 2 ^ (s.logratio + 1 - diploid_level) for s=segments if !isnan(s.depth))
			statuses = obs_cn_statuses
			if size(obs_cn_statuses, 1)==0
				statuses = hcat([[maj, tot] for tot=1:4 for maj=ceil(Int, tot/2):tot]...)'
			end
			for (maj_cn, tot_cn) in eachrow(statuses)
				l = expected_hsaf(cancer_fraction, maj_cn, tot_cn, 
							depth = round(Int, k*(cancer_fraction*tot_cn + (1-cancer_fraction)*2)))
				if l > ymax; continue; end
				if tot_cn > 4 && (sum(seg_tot_cns .== tot_cn) == 0 || 
								  maj_cn > maximum(seg_maj_cns[seg_tot_cns .== tot_cn])+1)
					continue
				end
				label = "$(tot_cn-maj_cn)+$maj_cn"
				line_plot([0, CHR_STARTS[23]], fill(l, 2), color=CN_PALETTE[tot_cn], line_width=0.5)
				text(CHR_STARTS[25]+3e7, l, label, h_align="left", size=5, color=CN_PALETTE[tot_cn])
			end

			# Plot segment medians
			for i in 1:length(segments)
				s = segments[i]
				if CHR_STARTS[findone(CHR_NAMES, s.chromosome)] in CHR_STARTS[23:25]; continue; end;
				start, stop = CHR_STARTS[findone(CHR_NAMES, s.chromosome)]+s.start, CHR_STARTS[findone(CHR_NAMES, s.chromosome)]+s.stop
				hsaf_mean, std = s.hsaf, s.hsaf_stdev
				if isnan(std); continue; end
				color=CN_PALETTE[seg_tot_cns[i]]
				line_plot([start, stop], fill(hsaf_mean, 2), color=color, line_width=2)
				if show_quantile_thresholds
					# Horizontal line showing the threshold quantiles of the segment
					@assert hsaf_q_thresh<0.5 "quantile thresholds must be between 0 and 5"
					lq, uq = quantile(Normal(hsaf_mean, std), [hsaf_q_thresh, 1-hsaf_q_thresh])
					line_plot(fill((start+stop)/2, 2), [lq, uq], color=color, line_width=1.2)
				end
			end

			# Scatter plot of hsaf values
			genome_scatter_plot(hetz_snp_fraction.chromosome,
								hetz_snp_fraction.position,
								hetz_snp_fraction.value)

			ylim(0.5, ymax)
			ylabel("Hetz SNP allele fraction")
			xlabel("Chromosome")

		# HSAF + CN scatter plot
		subplot((2,1), rows=2, cols=2)
			# Plot segment points
			x_seg = [logratio_to_cn(s.logratio - diploid_level, cancer_fraction, 2) for s=segments] # raw segment CN values
			y_seg = [s.hsaf for s=segments] # segment hsaf values
			s = [1+((s.stop-s.start) / 1e7) for s=segments] # point sizes scaled with segment lengths
			c = [CN_PALETTE[cn] for cn=seg_tot_cns]
			scatter_plot(x_seg, y_seg, color=c, size=s)

			# Plot expected levels as circles
			if !adaptive_axis_limits
				tots = 0:min(floor(Int, maximum(x_seg)), 14)
			elseif size(obs_cn_statuses, 1)!=0
				tots = 0:maximum(obs_cn_statuses)
			else
				tots = 0:4
			end
			x_exp = [tot for tot=tots for maj=ceil(Int, tot/2):tot]
			y_exp = [expected_hsaf(cancer_fraction, maj, tot, 
					 depth=round(Int, k*(cancer_fraction*tot + (1-cancer_fraction)*2)))
					 for tot=tots for maj=ceil(Int, tot/2):tot]
			c = [lighter_color(CN_PALETTE[cn]) for cn in x_exp]
			scatter_plot(x_exp, y_exp, border_color=c, fill_color="none", size=150)
			
			if adaptive_axis_limits && size(obs_cn_statuses, 1)!=0
				xmax = maximum(obs_cn_statuses[:,2])*1.1
			else
				xmax = maximum(x_seg[[s.stop-s.start>2e6 for s=segments]])*1.05
			end
			y_exp = y_exp[x_exp.<xmax,:]
			xlim(-0.2, xmax)
			ylim(0.5, max(maximum(y_exp), maximum(filter(!isnan, y_seg)))*1.05)
			
			ylabel("Hetz SNP allele fraction")
			xlabel("Copy number")
			grid(b=false)

		# MAF distributions
		subplot((2,2), rows=2, cols=2)
		maf_kdes(cancer_fraction, m_tot, m_alt, m_tot_cns; min_muts=min_visualized_muts)

		# title
		if title_prefix==nothing
			suptitle(@sprintf("Cancer fraction: %d%%, Diploid level: %.2f", cancer_fraction*100, diploid_level))
		else
			suptitle(@sprintf("%s - Cancer fraction: %d%%, Diploid level: %.2f", title_prefix, cancer_fraction*100, diploid_level))
		end
	end
end


# AVG MACN CLUSTERING
function assign_mutation_clusters(cancer_fraction, cluster_macn, cluster_size, 
						   alt_reads, total_reads, maj_cn, total_cn; min_mutations=1)
	@assert !any(all(total_cn.==0, dims=2)) "Mutation(s) with region total CN = 0 in all samples not allowed"
	@assert length(cancer_fraction) == size(cluster_macn, 2) == 
				size(alt_reads, 2) "Number of samples does not match"
	M = size(alt_reads, 1)
	S = length(cancer_fraction)
	C = size(cluster_macn, 1)
	probabilities = zeros(M, S, C)
	for m in 1:M
		m_alt = alt_reads[m, :]
		m_tot = total_reads[m, :]
		seg_maj = maj_cn[m, :]
		seg_tot = total_cn[m, :]
		cmb = binomial.(big.(m_tot), big.(m_alt))
		for c in 1:C
			# Check if cluster is possible for the mutation based on segment copy number status
			if any(cluster_macn[c, :] .> seg_maj)
				probabilities[m, :, c] .= 0
				continue
			end
			for s in 1:S
				maf_exp = expected_maf(cancer_fraction[s], cluster_macn[c, s], seg_tot[s], 2)
				p_obs = cmb[s] * maf_exp^m_alt[s] * (1-maf_exp)^(m_tot[s]-m_alt[s])
				probabilities[m, s, c] = p_obs
			end
		end
	end
	if cluster_size==nothing
		scores = prod(probabilities, dims=2)
		clusters = mapslices(argmax, scores, dims=3)[:]
		cluster_sizes = [sum(clusters.==c) for c=1:C]
		scores = mapslices(k -> k .* cluster_sizes, scores, dims=3)
	else
		probabilities .*= reshape(cluster_size, (1,1,length(cluster_size)))
		scores = prod(probabilities, dims=2)
	end
	clusters = mapslices(argmax, scores, dims=3)[:]
	clusters[mapslices(s -> all(s.==0), scores, dims=3)[:]] .= 0
	return clusters, scores[:,1,:]
end
function clust_macn_to_cf(macn; combine_dist=0.1, dont_combine=[])

	# Generally avg MACNs should be integers above one, but e.g. when
	# there is a minor clone with a very low fraction, its CN differences
	# do not affect coverage enough to not assign a common CN to them.
	# In such case, in a region where the major clone has higher CN than
	# the minor clone, the average MACN of a mutation may be little shy
	# of an integer greater than one. Therefore, we round the avg MACN
	# to get the real MACN, with which we divide the avg MACN to get
	# cancer fraction
	cf = [f>=1 ? f/round(f) : f for f=macn]

	dist_m = [c2<=c1 ? 1.0 : sqrt(sum((cf[c1,:].-cf[c2,:]).^2))
				for c1=1:size(cf,1), c2=1:size(cf,1)]
	to_combine = collect.(Tuple.(findall(dist_m .< combine_dist)))
	filter!.(c->!(c in dont_combine), to_combine)
	p1=1
	while p1<length(to_combine)
		c1=1
		while c1<=length(to_combine[p1])
			p2=p1+1
			while p2<=length(to_combine)
				if to_combine[p1][c1] in to_combine[p2]
					to_combine[p1] = unique(vcat(to_combine[p1],
										to_combine[p2]))
					deleteat!(to_combine, p2)
				end
				p2+=1
			end
			c1+=1
		end
		p1+=1
	end
	filter!(g -> any(macn[g,:] .> 1.0), to_combine)
	keep = setdiff(1:size(cf,1), vcat(to_combine...))
	original_clusters = [[[k] for k=keep]..., to_combine...]
	for g in to_combine
		cf[g[1], :] .= mean(cf[g,:], dims=1)[:]
		push!(keep, g[1])
	end
	final_order = sortperm(keep)
	cf = cf[keep[final_order], :]
	cf[findmax(sum(cf, dims=2)[:])[2],:] .= 1.0 # force truncal cluster to 100%
	return cf, original_clusters[final_order]
end
function simulate_ccf(total_reads::Matrix, total_cn::Matrix, macn::Matrix, ccf::Matrix, 
					cancer_frac::Vector; cfdna=falses(length(cancer_frac)), 
					sample_names=nothing, plot_macn=false, read_len=150, 
					maf_thresh=0.1, nreads_thresh=5, axis_limits=nothing)

	# rows are mutations and columns are samples in total_reads, total_cn, ccf & macn

	cfdna_frag_len = 170
	cfdna_reads_to_frags = cfdna_frag_len / (2*read_len)
	cfdna_two_reads_p = (2*read_len-cfdna_frag_len) / cfdna_frag_len

	sim_ccf = zeros(size(total_reads,1), size(total_reads,2));
	for k in 1:size(total_reads,1)
		p = ccf[k,:] .* expected_maf.(cancer_frac, macn[k,:], total_cn[k,:], 2)
		fragment_dist = Binomial.(round.(Int, total_reads[k,:].*[s ? cfdna_reads_to_frags : 1 for s=cfdna]), p)
		alt = fill(-Inf, size(total_reads,2))
		while all(alt./total_reads[k,:] .< maf_thresh) || all(alt .< nreads_thresh)
			alt_frags = rand.(fragment_dist)
			double_read_dist = Binomial.(alt_frags, cfdna .* cfdna_two_reads_p)
			alt = alt_frags .+ rand.(double_read_dist)
		end
		mccf = mutant_allele_cn.(alt./total_reads[k,:], cancer_frac, total_cn[k,:], 2)
		if !plot_macn
			mccf ./= macn[k,:]
		end
		sim_ccf[k,:] .= mccf
	end

	uniq_ccf = unique(ccf, dims=1);
	if size(uniq_ccf, 1) > length(GROUP_PALETTE)
		color = RGB(0)
	else
		clusters = Dict([ccf => k for (k, ccf) in enumerate(eachrow(uniq_ccf))]);
		cluster = [clusters[ccf[k,:]] for k in 1:size(total_reads,1)];
		color = GROUP_PALETTE[cluster]
	end

	nplots = binomial(size(total_reads,2), 2)
	figure(size=(5, nplots*5)) do
		p=0
		for s1 in 1:size(total_reads,2)-1, s2 in s1+1:size(total_reads,2)
			subplot(p+=1, nplots)
			scatter_plot(sim_ccf[:,s1], sim_ccf[:,s2], color=color, size=5)
			if axis_limits==nothing
				xlim(-0.05, sort(sim_ccf[:,s1])[end-ceil(Int, size(sim_ccf,1)/500)]+0.1)
				ylim(-0.05, sort(sim_ccf[:,s2])[end-ceil(Int, size(sim_ccf,1)/500)]+0.1)
			else
				xlim(axis_limits[p][1]...)
				ylim(axis_limits[p][2]...)
			end
			if sample_names!=nothing
				xlabel(sample_names[s1])
				ylabel(sample_names[s2])
			end
		end
	end
end
function mutation_genome_positions_plot(cancer_fraction::Real, diploid_level::Real, 
										coverage_logratio::GenomeTrack, hetz_snp_fraction::GenomeTrack,
										segments::Vector{Segment}, seg_tot_cns::Vector, seg_maj_cns::Vector,
										m_chr::Vector{String}, m_pos::Vector{T} where T <: Integer, 
										clusters::Vector{T} where T <: Integer; 
										tick_width=2e6, sample_name=nothing, omit_xy=false)							
	
	# Draw mutations as vertical ticks on top of genome tracks
	#   grouped and colored by clusters
	function draw_mut_positions()
		ylims = get_ylim()
		ticklength = (ylims[2] - ylims[1]) / 20
		new_ymax = ylims[2] + maximum(clusters) * ticklength
		rectangle(0, CHR_STARTS[25], ylims[2], new_ymax, color=RGB(255), zorder=5, edge_width=3)
		integer = seg_tot_cns.!=nothing
		for r in runs(integer)
			if integer[r[1]]; continue; end
			start = CHR_STARTS[findone(CHR_NAMES, segments[r[1]].chromosome)] + segments[r[1]].start
			stop = CHR_STARTS[findone(CHR_NAMES, segments[r[end]].chromosome)] + segments[r[end]].stop
			rectangle(start, stop, ylims[2], new_ymax, color=RGB(230), zorder=5)
		end
		if omit_xy
			rectangle(CHR_STARTS[23], CHR_STARTS[25], ylims[2], new_ymax, color=RGB(230), zorder=5)
		end

		for c in 1:maximum(clusters)
			top = new_ymax-(c-1)*ticklength
			bottom = new_ymax-c*ticklength
			muts = findall(clusters.==c)
			x = zeros(2 + 4*length(muts))
			y = zeros(2 + 4*length(muts))
			y[1]=bottom; y[end]=bottom
			for (i, m) in enumerate(muts)
				centre = CHR_STARTS[findone(CHR_NAMES, m_chr[m])] + m_pos[m]
				left = centre - tick_width/2
				right = centre + tick_width/2
				x[4*i-2]=left;  y[4*i-2]=bottom
				x[4*i-1]=left;  y[4*i-1]=top
				x[4*i  ]=right; y[4*i  ]=top
				x[4*i+1]=right; y[4*i+1]=bottom
			end
			area_plot(x, y, xbottom=bottom, zorder=6,
						color=(c <= length(GROUP_PALETTE) ? GROUP_PALETTE[c] : RGB(0)))
		end

		ylim(NaN, new_ymax)
	end

	# set max cn to max number in color palette
	max_cn = maximum(k for k=keys(CN_PALETTE) if typeof(k)==Int)
	for i in 1:length(seg_tot_cns)
		if seg_tot_cns[i]==nothing; continue; end
		if seg_tot_cns[i]>max_cn
			seg_tot_cns[i], seg_maj_cns[i] = nothing, nothing
		end
	end

	figure(size=(10,4)) do
		# Logratios
		subplot(1, 2)
		genome_scatter_plot(coverage_logratio.chromosome, 
									coverage_logratio.position, 
									coverage_logratio.value)
		genome_segments_plot(segments, [s.logratio for s=segments], 
			color=[CN_PALETTE[k] for k=seg_tot_cns], line_width=3)
		
		ymax, ymin = -Inf, Inf
		for cn in unique(seg_tot_cns)
			if cn==nothing; continue; end
			l = cn_to_logratio(cn, cancer_fraction, 2) + diploid_level
			if l<ymin; ymin=l-0.1; end
			if l>ymax; ymax=l+0.1; end
			line_plot([0, CHR_STARTS[25]], fill(l, 2), 
				color=CN_PALETTE[cn], line_width=0.5)
			text(CHR_STARTS[25]+3e7, l, cn, h_align="left", 
				size=7, color=CN_PALETTE[cn])
		end

		if ymax==-Inf
			ylim(minimum(coverage_logratio.values), maximum(coverage_logratio.values))
		else
			ylim(ymin, ymax)
		end
		draw_mut_positions()
		if sample_name!=nothing
			title(sample_name)
		end
		ylabel("Coverage logratio")
		xlabel("Chromosome")
		grid(false, axis="y")

		# HSAFs
		subplot(2, 2)
		genome_scatter_plot(String.(hetz_snp_fraction.chromosome), 
									hetz_snp_fraction.position, 
									hetz_snp_fraction.value)
		genome_segments_plot(segments, [s.hsaf for s=segments], 
			color=[CN_PALETTE[k] for k=seg_tot_cns], line_width=3)

		ymax = -Inf
		k = median(s.depth / 2^(s.logratio+1-diploid_level) for s=segments if !isnan(s.depth))
		for (maj_cn, tot_cn) in unique(zip(seg_maj_cns, seg_tot_cns))
			if maj_cn==nothing; continue; end
			l = expected_hsaf(cancer_fraction, maj_cn, tot_cn, 
					depth=round(Int, k*(cancer_fraction*tot_cn + (1-cancer_fraction)*2)))
			if l>ymax; ymax=l+0.1; end
			line_plot([0, CHR_STARTS[25]], fill(l, 2), color=CN_PALETTE[tot_cn], line_width=0.5)
			text(CHR_STARTS[25]+3e7, l, "$maj_cn+$(tot_cn-maj_cn)", h_align="left", 
				size=7, color=CN_PALETTE[tot_cn])
		end

		if ymax==-Inf
			ylim(NaN, maximum(hetz_snp_fraction.values))
		else
			ylim(NaN, ymax)
		end
		draw_mut_positions()
		grid(false, axis="y")
		ylabel("Hetz SNP allele fraction")
		xlabel("Chromosome")
	end
end
# dirichlet



# Since in each iteration of the clustering process mutations are stochastically
# assigned to clusters, many clusters get no mutations assigned to them (i.e.
# are eliminated). If there are N clusters left after an iteration, this
# function relabels all remaining clusters to have identifiers 1...N.
function compactify_clusters!(cluster::Vector{Int})
	C = 0
	relabel = zeros(Int, maximum(cluster))
	for m in 1:length(cluster)
		c = cluster[m]
		if relabel[c] > 0; continue; end
		relabel[c] = C + 1   # Rename cluster c -> C + 1
		C += 1
	end
	cluster .= relabel[cluster]
	return C    # Return the number of clusters after compactification
end
function bayesian_clustering(cancer_fraction, alt_reads, total_reads, total_cn;
	iterations=3000, initial_clusters=1000, dispersion_factor=0.25,
	silent=false)
	
	#iterations=discovery_iterations+merge_iterations
	S = length(cancer_fraction)
	M = size(alt_reads, 1)
	#maf = alt_reads ./ total_reads
	macn = fill(NaN, M, S)
	for s in 1:S
		macn[:, s] = mutant_allele_cn.(alt_reads[:, s] ./ total_reads[:, s],
			cancer_fraction[s], total_cn[:, s], 2)
	end

	# To start, we create a user-defined number of initial clusters, each placed
	# on a randomly selected mutation. All mutations are assigned to their
	# nearest initial cluster.
	centroids = macn[randperm(M)[1:initial_clusters], :]
	cluster = map(1:M) do m
		argmin([sum((macn[m, :] .- centroids[c, :]).^2)
			for c in 1:initial_clusters])
	end

	for iter in 1:iterations

		# At the start of every new iteration, we remove empty clusters and
		# compactify cluster assignments such that all clusters 1..C contain
		# at least one mutation
		C = compactify_clusters!(cluster)

		# Calculate the number of mutations in each cluster
		num_cluster_muts = StatsBase.counts(cluster, C)

		# Calculate cluster centroid coordinates (in MACN space)
		centroids = fill(NaN, C, S)
		for s in 1:S
			macn_sums = zeros(C)
			for (m, c) in enumerate(cluster)
				macn_sums[c] += macn[m, s]
			end
			centroids[:, s] .= macn_sums ./ num_cluster_muts
		end

		if !silent && (iter % 50) == 0
			@printf(stderr, "Iteration %d (%d clusters)...\n", iter, C)
		end

		# Randomly assign each mutation to one of the existing clusters, with
		# the assignment probabilities determined based on the relative
		# likelihoods of the mutation arising from each cluster. Or assign
		# the mutation into its own new cluster, at probability alpha.
		for m in 1:M
			# Calculate likelihood of mutation belonging to each MACN cluster
			cluster_lh = ones(C)
			for c in 1:C
				for s in 1:S
					mut_cn = centroids[c, s]
					if mut_cn > 1.1 * total_cn[m, s]
						cluster_lh[c] = 0 # MACN cannot be higher than region CN
						break
					end
					mut_cn = min(mut_cn, total_cn[m, s])

					emaf = expected_maf(cancer_fraction[s], mut_cn, total_cn[m, s], 2)
					dist = Binomial(round(Int, dispersion_factor * total_reads[m, s]), emaf)
					cluster_lh[c] *= pdf(dist, round(Int, dispersion_factor * alt_reads[m, s]))
				end
			end

			# Calculate assignment probabilities (likelihood * cluster size)
			weighted_lh = cluster_lh .* num_cluster_muts

			# Sample from a categorical distribution to determine the new
			# cluster assignment
			total = sum(weighted_lh)
			if total == 0
				# If all clusters have a likelihood of zero (it can happen),
				# give all clusters an identical assignment probability.
				assignment_p = fill(1 / C, C)
			else
				assignment_p = weighted_lh ./ total
			end

			cluster[m] = rand(Categorical(assignment_p))
		end
	end

	C = compactify_clusters!(cluster)
	println("Finished with $C clusters.")

	# Calculate cluster centroids
	cluster_macn = fill(NaN, C, S)
	for c in 1:C
		cluster_muts = findall(cluster .== c)
		for s in 1:S
			cluster_macn[c, s] = mean(macn[cluster_muts, s])
		end
	end

	# force cluster closest to 1 to 1
	cluster_macn[argmin(abs.(sum(cluster_macn .- 1, dims=2)[:])), :] .= 1.0

	# order clusters from highest to lowest sum over macns of all samples
	cluster_order = sortperm(sum(cluster_macn, dims=2)[:], rev=true)
	cluster_macn = cluster_macn[cluster_order,:]
	cluster = sortperm(cluster_order)[cluster]

	return cluster_macn, cluster
end


# PHYLOGENETICS
function permuted_partitions_with_empties(s::AbstractVector, m::Int)
	ret = Vector{Vector{Int}}[]
	for n_empty in 0:m-1
		for non_empties in combinations(1:m, m-n_empty)
			for groups in partitions(s, m-n_empty)
				for order in permutations(1:m-n_empty)
					parts = fill(Int[], m)
					parts[non_empties] = groups[order]
					push!(ret, parts)
				end
			end
		end
	end
	return ret
end
function solve_phylogeny_matrix(fracs; max_fracs=ones(size(fracs, 2)))
	solutions = []
	C = size(fracs, 1)
	if C==0
		return [falses(0, 0)]
	end
	parent_combs = collect(combinations(1:C))
	filter!(p -> all(sum(fracs[p, :], dims=1)[:] .<= max_fracs), parent_combs)
	for parents in parent_combs
		children = setdiff(1:C, parents)
		Np = length(parents)
		Nc = length(children)
		if Nc==0
			push!(solutions, falses(C, C))
			continue
		end
		for sibling_groups in permuted_partitions_with_empties(children, Np)
			branch_solutions = [solve_phylogeny_matrix(fracs[c, :], max_fracs=fracs[p, :])
								for (p, c) in zip(parents, sibling_groups)]
			if any(length.(branch_solutions) .== 0); continue; end
			for branch_solution_comb in Base.product(branch_solutions...)
				solution = falses(C, C)
				for (p, c, branch_solution) in zip(parents, sibling_groups, branch_solution_comb)
					solution[c, p] .= true
					solution[c, c] = branch_solution
				end
				push!(solutions, solution)
			end
		end
	end
	return solutions
end
function phylo_matrix_to_chain(M)
	parents = zeros(Int, size(M,1))
	height = sum(M, dims=1)[:]
	for m in 1:size(M,1)
		papas = findall(M[m,:])
		if length(papas)==0; continue; end
		parents[m] = papas[findmin(height[papas])[2]]
	end
	return parents
end
function solve_phylogeny(fracs)
	phylogeny_matrices = solve_phylogeny_matrix(fracs)
	ret = []
	for parents in phylo_matrix_to_chain.(phylogeny_matrices)
		subclone_fracs = zeros(size(fracs))
		for m in 1:size(fracs, 1),  s in 1:size(fracs, 2)
			subclone_fracs[m, s] = fracs[m, s] - sum(fracs[findall(parents.==m), s])
		end
		push!(ret, (parents=parents, fractions=subclone_fracs))
	end
	return ret
end
function phylotree_plot(parents; sizes=nothing, x=0, y=0, width=1, height=nothing, line_width=10,
						h_align="left", v_align="bottom", color=RGB(0), dir=1)
	descendants(node, children) = vcat(children[node]..., 
		[descendants(child_node, children) for child_node=children[node]]...)
	branch_depth(node, parents; depth=1) = parents[node]==0 ? depth : 
		branch_depth(parents[node], parents, depth=depth+1)
	
	children = [findall(parents.==node) for node=1:length(parents)]
	descs = [descendants(node, children) for node=1:length(parents)]

	xcoords = zeros(length(parents))
	for node in 1:length(parents)
		middle_child_i = length(children[node])/2+0.5
		for (i, child_node) in enumerate(children[node])
			branch_width = sum([length(children[n])==0 ? 0 : length(children[n])-1
								for n=vcat(child_node, descs[child_node])])
			for child_node_right in children[node][i+1:end]
				xcoords[child_node_right] += branch_width+1
				for child_node_right_desc in descs[child_node_right]
					xcoords[child_node_right_desc] += branch_width+1
				end
			end
			shift = 0.0
			if i < middle_child_i
				shift += branch_width+0.5
			elseif i == middle_child_i && length(children[node])!=1
				shift += branch_width/2+0.5
			end
			nd = copy(node)
			while true
				xcoords[nd] += shift
				nd = parents[nd]
				if nd==0 || length(children[nd]) != 1; break; end
			end
		end
	end

	if sizes!=nothing
		depths = copy(sizes)
		for node in 1:length(parents), desc_node in descs[node]
			depths[desc_node] += sizes[node]
		end
	else
		sizes = ones(length(parents))
		depths = [branch_depth(node, parents) for node=1:length(parents)]
	end

	depths = -depths .+ maximum(depths)
	if height==nothing
		height = (maximum(depths) .+ sizes[findone(parents, 0)])
		height_multiplier = 1
	else
		height_multiplier = 1 / (maximum(depths) .+ sizes[findone(parents, 0)]) .* height
	end
	sizes = Float64.(sizes) .* height_multiplier
	depths = depths .* height_multiplier
	if !all(xcoords.==0)
		xcoords = xcoords ./ maximum(xcoords) .* width
	end
	if (dir==1 && h_align=="right") || (dir==2 && v_align=="top")
		xcoords .-= maximum(xcoords)
	elseif (dir==1 && h_align=="centre") || (dir==2 && v_align=="centre")
		xcoords .-= maximum(xcoords)/2
	end
	if (dir==1 && v_align=="top") || (dir==2 && h_align=="right")
		depths .-= maximum(depths) + sizes[findone(parents, 0)]
	elseif (dir==1 && v_align=="centre") || (dir==2 && h_align=="centre")
		depths .-= (maximum(depths) + sizes[findone(parents, 0)])/2
	end

	xcoords .+= dir==1 ? x : y
	depths .+= dir==1 ? y : x

	figure() do
		for node in (1:length(parents))[sortperm(depths)]
			c = typeof(color)==RGB ? color : color[node]
			xy = fill(xcoords[node], 2), [depths[node]+sizes[node], depths[node]]
			line_plot((dir==1 ? xy : (height .- xy[2], xy[1]))..., 
					color=c, line_width=line_width, capstyle="round")
			if length(children[node]) > 1
				xy = xcoords[children[node]][[1, end]], fill(depths[node], 2)
				line_plot((dir==1 ? xy : (height .- xy[2], xy[1]))..., 
						color=c, line_width=line_width, capstyle="round")
			end
		end
	end
end


# SUBCLONE COPY NUMBER DECONVOLUTION

# used to give compiler better information on data types
@generated function ntuple_type(::Val{N}, ::Val{T}) where {N,T}
	:(NTuple{$N, $T})
end
function cn_deconvolution(cancer_fraction::AbstractVector, diploid_level::AbstractVector,
						   segments::Matrix{Segment}, sc_fracs::Matrix{Float64}, 
						   sc_parents::AbstractVector; lr_const=1, hsaf_const=0.15, 
						   cna_const=0.05, max_cn=7, max_0_seg_length=2e7, 
						   n_clones_as_type::Val{N_CLONES}=Val(size(sc_fracs, 1))) where N_CLONES # <- used by the above defined ntuple_type() to give compiler more information user should not touch this!
	# Rules for model fit
	# 1. Minimize abs(expected coverage - observed coverage)
	# 2. Minimize abs(expected hsaf - observed hsaf)
	# 3. Minimize number of CNAs needed for CN statuses
	# 4. No large areas with 0 copies allowed
	# 5. No allele specific copy number can be gained from 0

	@assert length(cancer_fraction) == length(diploid_level) == 
										size(segments,2) == size(sc_fracs, 2)
	@assert size(sc_fracs, 1) == length(sc_parents)
	@assert all(getfield(segments[seg, 1], field) == getfield(segments[seg, sample], field) 
				for sample=2:size(segments,2) for seg=1:size(segments,1)
				for field=[:chromosome, :start, :stop])
	@assert all(sum(sc_fracs, dims=1).-1 .< 0.00001)
	n_clones = size(sc_fracs, 1);
	n_segs = size(segments, 1);

	if deviated_hsaf_table==nothing
		# run this once just to load the table
		expected_hsaf(1.0, 1, 2, depth=5, use_table=true)
	end
	hsaf_dev_table::Matrix{Float64} = copy(deviated_hsaf_table)

	# define product iterator type for compiler
	T = ntuple_type(Val(N_CLONES), Val(UnitRange{Int64}))
	@assert length(T.parameters) == N_CLONES
	# clone CN state search space
	tot_space = Base.product(fill(0:max_cn, n_clones)...)::Base.Iterators.ProductIterator{T};
	n_alts = [*((tots.+1)...) for tots=tot_space];
	cns_space = fill(-1, sum(n_alts), n_clones, 2);
	minimum_cnas = zeros(Int, sum(n_alts))
	r=0;
	for tot in tot_space, alt in Base.product([0:c for c=tot]...)::Base.Iterators.ProductIterator{T}
		r+=1; valid=true
		for cl in 2:n_clones
			t, m = tot[cl], alt[cl]
			pt, pm = tot[sc_parents[cl]], alt[sc_parents[cl]]
			if (pt==0 && t!=0) || (pm==0 && m!=0) || (pt==pm && t!=m)
				valid=false; break
			end
			minimum_cnas[r] += abs(pm-m) + abs(pt-pm-t+m)
		end
		if !valid; continue; end
		cns_space[r, :, 1] .= tot
		cns_space[r, :, 2] .= alt
	end
	valid = cns_space[:,1,1].!=-1
	minimum_cnas = minimum_cnas[valid]
	cns_space = cns_space[valid,:,:]
	order = sortperm(minimum_cnas)
	minimum_cnas = minimum_cnas[order]
	cns_space = cns_space[order,:,:]

	# calculate expected coverage log ratios and hsafs for cn states
	tots_exp =  sum(sc_fracs .* permutedims(cns_space[:,:,[1]], (2, 3, 1)), dims=1)[1,:,:]';
	majs_exp = sum(sc_fracs .* permutedims(cns_space[:,:,[2]], (2, 3, 1)), dims=1)[1,:,:]';
	logr_exp = [cn_to_logratio(tots_exp[c, s], cancer_fraction[s], 2) + diploid_level[s]
				for c=1:size(tots_exp, 1), s=1:size(tots_exp, 2)];
	hsaf_exp = [0.5 + abs(0.5 - expected_hsaf(cancer_fraction[s], majs_exp[c,s], tots_exp[c,s]))
				for c=1:size(tots_exp, 1), s=1:size(tots_exp, 2)];
	hsaf_exp_corr_i = [round(Int, (f-0.5)*1000)+1 for f=hsaf_exp]

	# initialize result arrays
	tot_cns = Array{Union{Nothing, Int}}(nothing, n_segs, n_clones)
	alt_cns = Array{Union{Nothing, Int}}(nothing, n_segs, n_clones)
	errs = fill(Inf, n_segs);

	# find best cn combination for each segment separately
	for seg_i in 1:n_segs

		segs = segments[seg_i,:]
		if any(isnan(s.depth) || isnan(s.hsaf_stdev) for s=segs)
			continue
		end

		# observed segment data
		is_short = segs[1].stop - segs[1].start > max_0_seg_length;
		logr_obs = [s.logratio for s=segs];
		hsaf_obs = [s.hsaf for s=segs];
		hsaf_obs_std = [s.hsaf_stdev for s=segs];
		depth_corr_i = min.(round.(Int, s.depth for s=segs), size(hsaf_dev_table, 1));

		# keep track of best search result
		best_err = Inf;
		best_cns = nothing;
		best_i=0

		# exhaustive search
		for i in 1:size(cns_space, 1)
			err = cna_const * minimum_cnas[i];
			if err > best_err; break; end;
			err += lr_const * maximum(abs(exp - obs) for (exp,obs)=zip(logr_exp[i,:], logr_obs));
			if err > best_err; continue; end;
			err += hsaf_const * maximum(abs(obs - hsaf_dev_table[di, hi]) / std
				for (obs, std, hi, di)=zip(hsaf_obs, hsaf_obs_std, hsaf_exp_corr_i[i,:], depth_corr_i));
			if err > best_err; continue; end;
			if is_short && any(t==0 for t=cns_space[i,:,1]); continue; end;
			best_err = err;
			best_cns = cns_space[i,:,:];
			best_i=i
		end

		errs[seg_i] = best_err;
		tot_cns[seg_i, :] .= best_cns[:,1];
		alt_cns[seg_i, :] .= best_cns[:,2];
	end

	return tot_cns, alt_cns, errs
end
function cn_deconvolution_plot(cancer_fraction::AbstractVector, 
							   diploid_level::AbstractVector,
							   coverage_logratio::Vector{GenomeTrack{T}} where T <: Real, 
							   hetz_snp_fraction::Vector{GenomeTrack{T}} where T <: Real,
							   segments::Matrix{Segment}, sc_fracs::AbstractMatrix, 
							   sc_parents::Vector, tot_cns::AbstractArray,
							   alt_cns::AbstractArray;
							   bad_deconvolutions=falses(size(tot_cns, 1)),
							   plot_sex_chrs=false, sample_names=nothing, clone_colors=nothing)
	@assert length(cancer_fraction) == length(diploid_level) == 
			length(coverage_logratio) == length(hetz_snp_fraction)
	@assert all(sum(sc_fracs, dims=1).-1 .< 0.00001)

	if sample_names != nothing
		@assert length(sample_names) == length(cancer_fraction)
	end

	if clone_colors==nothing
		clone_colors = GROUP_PALETTE[1:length(sc_parents)]
	end

	n_samples = length(cancer_fraction)
	# n_segs = size(segments,1)
	n_clones = size(sc_fracs, 1);
	last_chr_i = plot_sex_chrs ? 24 : 22

	figure(size=(15, 3.5*(n_clones+n_samples))) do
		gridrows = 4*n_clones + 8*n_samples + 1
		clone_rspan = 4
		sample_rspan = 4
		clone_marker_rspan = 4
		subclone_fractions_rspan = 8
		gridcols = 30
		genome_colspan = 29
		r = 1
		
		# Subclone allele-specific copy number profiles
		for sc_i in 1:n_clones
			# Raw copy numbers
			subplot((r, 1), rows=gridrows, cols=gridcols, rowspan=clone_rspan, colspan=genome_colspan)

				for (seg_i, s) in enumerate(segments[:,1])
					if isnothing(tot_cns[seg_i, 1]); continue; end

					tot_cn = tot_cns[seg_i, sc_i]
					alt_cn = alt_cns[seg_i, sc_i]
					color = CN_PALETTE[tot_cn]

					# draw segment
					s_chr_start = CHR_STARTS[findone(CHR_NAMES, s.chromosome)]
					
					area_plot([s_chr_start+s.start, s_chr_start+s.stop], [tot_cn+0.2, tot_cn+0.2], ybottom=tot_cn)
					area_plot([s_chr_start+s.start, s_chr_start+s.stop], [alt_cn, alt_cn], color=RGB(76,204,255))
					area_plot([s_chr_start+s.start, s_chr_start+s.stop], [tot_cn, tot_cn], ybottom=alt_cn, color=RGB(243,128,171))

					# mark bad deconvolutions with "x"
					if bad_deconvolutions[seg_i]
						middle = ((s_chr_start+s.start)+(s_chr_start+s.stop))/2
						@suppress_err scatter_plot([middle], [tot_cn], size=15, marker="x", color=CN_PALETTE["bad"], zorder=10)
					end
					
				end

				genome_plot_config(CHR_NAMES[1:last_chr_i], CHR_SIZES[1:last_chr_i])

				PyPlot.gca().set_xticklabels([], minor=true)
				xlim(0, CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i])
				ylabel("Copy number")
			
			# Subclone marker
			subplot((r, gridcols), rows=gridrows, cols=gridcols, rowspan=clone_marker_rspan);
				rectangle(0, 1, 0, 1, color=clone_colors[sc_i])
				xticks([], spine=false)
				yticks([], spine=false)

			r+=clone_rspan

		end

		# Sample LR scatters and expected segment levels
		for si in 1:n_samples
			# Coverage
			subplot((r, 1), rows=gridrows, cols=gridcols, rowspan=sample_rspan, colspan=genome_colspan)
				# coverage logratio scatterplot
				genome_scatter_plot(string.(coverage_logratio[si].chromosome),
											coverage_logratio[si].position,
											coverage_logratio[si].value)
				for (seg_i, s) in enumerate(segments[:,si])
					if isnothing(tot_cns[seg_i, 1]); continue; end

					# if there are no CNAs between subclones, set segment color according to cn
					if all((tot_cns[seg_i, :] .== tot_cns[seg_i, 1]) &
						   (alt_cns[seg_i, :] .== alt_cns[seg_i, 1]))
						cn = tot_cns[seg_i, 1]
						color=CN_PALETTE[cn]
					# otherwise set color to green
					else
						cn = sum(sc_fracs[:,si] .* tot_cns[seg_i, :])
						color = CN_PALETTE["subclonal"]
					end

					# draw expected segment level
					s_chr_start = CHR_STARTS[findone(CHR_NAMES, s.chromosome)]
					clr = cn_to_logratio(cn, cancer_fraction[si], 2) + diploid_level[si]
					line_plot([s_chr_start+s.start, s_chr_start+s.stop], fill(clr, 2), 
								color=color, line_width=2)

					# mark bad deconvolutions with "x"
					if bad_deconvolutions[seg_i]
						middle = ((s_chr_start+s.start)+(s_chr_start+s.stop))/2
						@suppress_err scatter_plot([middle], [clr], size=15, marker="x", zorder=10, color=CN_PALETTE["bad"])
					end
				end

				# define ymin & ymax
				drawn_cns = filter(cn -> cn!=nothing, tot_cns[:])
				min_drawn, max_drawn = minimum(drawn_cns), maximum(drawn_cns)
				min_obs_lvl = cn_to_logratio(min_drawn, cancer_fraction[si], 2)+diploid_level[si]
				min_ref_lvl = cn_to_logratio(min_drawn==0 ? 1 : min_drawn-1, cancer_fraction[si], 2)+diploid_level[si]
				max_ref_lvl = cn_to_logratio(max_drawn+2, cancer_fraction[si], 2)+diploid_level[si]
				sex_chr_filter = plot_sex_chrs ? trues(length(coverage_logratio)) : !occursin.(r"^.+(X|Y)", coverage_logratio[si].chromosome)
				min_lr_point = minimum(coverage_logratio[si].value[sex_chr_filter])
				max_lr_point = maximum(coverage_logratio[si].value[sex_chr_filter])
				ymin = max((min_obs_lvl+min_ref_lvl)/2, min_lr_point)
				ymax = min(max_ref_lvl, max_lr_point)

				ylim(ymin, ymax)
				xlim(0, CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i])
				PyPlot.gca().set_xticklabels([], minor=true)
				ylabel("Coverage logratio")
				grid(false, axis="y")

			# Subclone marker
			subplot((r, gridcols), rows=gridrows, cols=gridcols, rowspan=subclone_fractions_rspan);
				stacked_bar_plot(sc_fracs[end:-1:1,[si]]', colors=clone_colors[end:-1:1])
				xlim(0.8, 1)
				text(1.12, 0.5, sample_names==nothing ? "Sample-$si" : sample_names[si], rotation=-90)
				xticks([], spine=false)
				yticks([], spine=false)

			r+=sample_rspan

			# HSAF
			subplot((r, 1), rows=gridrows, cols=gridcols, rowspan=sample_rspan, colspan=genome_colspan)
				# hsaf scatter plot
				genome_scatter_plot(string.(hetz_snp_fraction[si].chromosome),
											hetz_snp_fraction[si].position,
											hetz_snp_fraction[si].value)
				
											
				for (seg_i, s) in enumerate(segments[:,si])
					if isnothing(tot_cns[seg_i, 1]); continue; end

					# find expected levels
					if all((tot_cns[seg_i, :] .== tot_cns[seg_i, 1]) &
						   (alt_cns[seg_i, :] .== alt_cns[seg_i, 1])) # if there are no CNAs between subclones
						cn_tot = tot_cns[seg_i, 1]
						cn_alt = alt_cns[seg_i, 1]
						lvl = expected_hsaf(cancer_fraction[si], cn_alt, cn_tot, 
								depth=s.depth, iterations=10_000)
						# set color according to cn
						color=CN_PALETTE[cn_tot]
					else
						cn_tot = tot_cns[seg_i, :]
						cn_alt = alt_cns[seg_i, :]
						cn_tot_avg = sum(sc_fracs[:,si] .* cn_tot);
						cn_alt_avg = sum(sc_fracs[:,si] .* cn_alt);
						lvl = expected_hsaf(cancer_fraction[si], cn_alt_avg, cn_tot_avg,
									depth=s.depth, iterations=10_000);
						# set color to green
						color = CN_PALETTE["subclonal"]
					end
					if lvl < 0.5
						lvl = 1 - lvl
					end

					# draw expected segment level
					s_chr_start = CHR_STARTS[findone(CHR_NAMES, s.chromosome)]
					line_plot([s_chr_start+s.start, s_chr_start+s.stop], fill(lvl, 2), 
								color=color, line_width=2)
					
					# mark bad deconvolutions with "x"
					if bad_deconvolutions[seg_i]
						middle = ((s_chr_start+s.start)+(s_chr_start+s.stop))/2
						@suppress_err scatter_plot([middle], [lvl], size=15, marker="x", color=CN_PALETTE["bad"], zorder=10)
					end

				end

				ylabel("HSAF")
				xlim(0, CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i])
				ylim(0.5, NaN)
				PyPlot.gca().set_xticklabels([], minor=true)
				grid(false, axis="y")

				r+=sample_rspan
		end

		# chromosome labels
		subplot((r, 1), rows=gridrows, cols=gridcols, colspan=genome_colspan)
			chr_labels = [r"(19|21|M)$" in chr ? "" : replace(chr, "chr", "")
						for chr in CHR_NAMES[1:last_chr_i]]
			chr_middles = cumsum(CHR_SIZES) .- CHR_SIZES ./ 2
			for (txt, pos) in zip(chr_labels, chr_middles)
				text(pos, 1, txt, v_align="top")
			end
			xlim(0, CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i])
			xticks([], spine=false)
			yticks([], spine=false)
			xlabel("Chromosome")
	end
end


# CLONAL EXPANSION + TUMOR BURDEN PLOT
function fit_exponential_function(x::AbstractVector, y::AbstractVector)
	@assert length(x) == length(y) == 2
	@assert !any(y.==0) "Exponential function can not go to 0"
	c = log(y[2]/y[1]) / (x[2]-x[1])
	a = y[1] / exp(c*x[1])
	f(x) = a*exp(c*x)
	return f
end
function branch_depth(node; children)
	if length(children[node])==0
		return 1
	else
		return 1+maximum(branch_depth.(children[node], children=children))
	end
end
function clonal_expansion_plot(psa, pt, scf, st, parents; zerostate_start=-50, smooth=false)
	@assert length(psa) == length(pt)
	@assert issorted(pt)
	
	# all parents must be before their children in "parents", i.e., it must be sorted
	children = [findall(parents.==sc) for sc in 1:length(parents)]
	@assert all([all(children[p] .> p) for p in 1:length(children)])

	# create 'beginning' timepoint and between-sample timepoints for visualizing clone births and deaths
	st2 = Float64.(hcat(st, st)'[:])
	st2 = vcat(st2[1]+zerostate_start, st2)
	for interval in 2:2:length(st2)-1
		st2[interval] = (st2[interval-1]+st2[interval+1])/2
	end
	st2 = vcat(fill(st2', length(parents))...)

	# create subclone fractions for beginning and between-sample timepoints using exponential interpolation
	scf2 = scf[:, hcat(1:size(scf, 2), 1:size(scf, 2))'[:]]
	scf2 = hcat(vcat(1.0, fill(0.0, size(scf,1)-1)), scf2)
	for c in 1:size(scf2, 1), i in 2:2:size(scf2, 2)-1
		if scf2[c, i-1]==0 && scf2[c, i+1]==0
			continue
		elseif scf2[c, i-1]==0 || scf2[c, i+1]==0
			# real exponential interpolation is not possible becaus the function never goes to 0
			#   we set the start of the function to 0.15 times adjacent level
			#   then we strech the line end to zero
			f = scf2[c, [i-1, i+1]]
			f[f.==0] .= f[f.!=0] .* 0.15
			itp = fit_exponential_function([0,1], f).([0, 0.5, 1.0])
			endpoint_shift = scf2[c, i-1]==0 ? (itp[1], 0) : (0, itp[end]) 
			itp .-= range(endpoint_shift..., length=length(itp))
			scf2[c, i] = itp[2]
		else
			scf2[c, i] = fit_exponential_function([0,1], scf2[c, [i-1, i+1]])(0.5)
		end
	end

	# calculate lineage fractions
	linf = copy(scf2)
	for sc in (1:size(linf, 1))[sortperm(parents, rev=true)]
		if parents[sc]==0; continue; end
		linf[parents[sc],:] .+= linf[sc,:]
	end

	# evenly time clones' between-sample deaths and births
	#   all deaths occur halfway between samples or after the last new clone is born
	#   nested births are scatterd evenly between samples
	for tp in 1:size(scf, 2)
		tp1 = 2*tp-1
		tp2 = 2*tp+1
		interval = 2*tp

		# check if any clones die
		tp1_0 = scf2[:, tp1] .== 0 #note scf
		after0 = [all(r .== 0) for r in eachrow(scf2[:, tp2:end])]
		dies = findall(!tp1_0 & after0)
		scf2[dies, interval] .= 0

		# check if any clones are born
		tp2_0 = linf[:, tp2] .== 0
		before0 = [all(r .== 0) for r in eachrow(linf[:, 1:tp1])]
		is_born = findall(before0 & !tp2_0)
		scf2[is_born, interval] .= 0

		# time clone births evenly according to phylogeny
		for sc in is_born
			depth = branch_depth(sc, children=children) + 1
			if length(dies) != 0
				depth += 1
			end
			if parents[sc] in is_born
				base_time = st2[parents[sc], interval]
			else
				base_time = st2[1, tp1]
			end
			st2[sc, interval] = base_time + (st2[1, tp2] - base_time) / depth
		end

		# time clone deaths according to max depth of born clones
		if length(is_born) != 0
			depth = branch_depth(is_born[1], children=children) + 2
			st2[dies, interval] .= st2[1, tp2] - (st2[1, tp2] - st2[1, tp1]) / depth
		end
	end
	linf = nothing # not correct anymore after timings were changed (and not used)

	# interpolate psa
	time_grid = st[1]+zerostate_start:maximum(st2)
	psa = vcat(0, psa)
	pt = vcat(st[1]+zerostate_start, pt)
	if smooth
		psai = interpolate(pt, psa, SteffenMonotonicInterpolation()).(time_grid)
	else
		psai = LinearInterpolation(pt, psa)(time_grid)
	end

	# Interpolate subclone fractions exponentially
	# scfi = vcat([LinearInterpolation(st2[sc,:], scf2[sc,:])(time_grid) for sc in 1:size(scf2, 1)]'...)
	scfi = vcat(map(1:size(scf2, 1)) do sc
		interp = Float64[]
		tt = st2[sc,:]
		for tp in 2:length(tt)
			t = [tt % 1 == 0.0 ? tt-0.01 : tt for tt in tt[tp-1:tp]]
			if tp==length(tt); t[2]+=0.01; end
			t_grid = ceil(t[1]):floor(t[2])
			f = scf2[sc,tp-1:tp]
			if all(f .== 0)
				itp = LinearInterpolation(t, f)(t_grid)
				#todo? itp = fill(0, length(grid))
			elseif any(f .== 0)
				# real exponential interpolation is not possible becaus the function never goes to 0
				#   we set the start of the function to 0.15 times adjacent level
				#   then we strech the line end to zero
				ff = copy(f)
				ff[ff.==0] .= ff[ff.!=0] .* 0.15
				itp = fit_exponential_function(t, ff).(t_grid)
				endpoint_shift = f[1]==0 ? (itp[1], 0) : (0, itp[end]) 
				itp .-= range(endpoint_shift..., length=length(itp))
			else
				itp = fit_exponential_function(t, f).(t_grid)
			end
			interp = vcat(interp, itp)
		end
		return interp
	end'...)

	# calculate interpolated lineage fractions
	linfi = copy(scfi)
	for sc in (1:size(linfi, 1))[sortperm(parents, rev=true)]
		if parents[sc]==0; continue; end
		linfi[parents[sc],:] .+= linfi[sc,:]
	end
	c = 1 ./ linfi[1,:]
	linfi = mapslices(r -> r .* c, linfi, dims=2)
	scfi = mapslices(r -> r .* c, scfi, dims=2)
	
	psi_sc = mapslices(r -> r .* psai, linfi, dims=2)
	
	y_bounds = cat(zeros(size(psi_sc)...), vcat(fill(psai', size(psi_sc, 1))...), dims=3)
	first_decendant_drawn = falses(size(psi_sc, 1))
	figure() do
		for sc in 1:size(psi_sc, 1)
			parent = parents[sc]
			if parent == 0
				ylb = zeros(size(psi_sc, 2))
				yub = psai
			else
				ylb = y_bounds[parent,:,1]
				yub = ylb .+ psi_sc[sc, :]
				# Raise the lower bound temporarily when a clone is born to indicate whose its parent
				if !first_decendant_drawn[parent]
					start = searchsortedlast(linfi[sc,:], 0)
					
					sc_death = findlast(linfi[sc,start+1:end] .== 0)
					sc_death = sc_death==nothing ? size(linfi, 2) : sc_death
					next_sample = findone(time_grid, st[findfirst(st .> time_grid[start])])
					stop = Int(round(4/5 * next_sample)) # "temporarily" is 4/5 towards next sample

					parent_lin_psa = y_bounds[parent,start,2] - y_bounds[parent,start,1]
					next_sample_parent_lin_psa = y_bounds[parent,next_sample,2] - y_bounds[parent,next_sample,1]
					next_sample_sibling_frac = 1 - scfi[parent, next_sample] / linfi[parent, next_sample]
					max_raise = parent_lin_psa * next_sample_sibling_frac  / 2
					raise = nothing
					while true
						# raise_len = stop-start+1
						# raise = fit_exponential_function([1, raise_len], [1.0, 1/raise_len]).(1:raise_len)
						# raise = strech_line_end_to_zero(raise, mode="end")
						# raise .*= max_raise
						raise_len = stop-start
						raise = max_raise .* (raise_len:-1:0) ./ raise_len
						# stop the raise earlier if it would be more than the parent's fraction
						max_raise_i = scfi[parent, start:stop] .* psai[start:stop]
						if findfirst(raise .> max_raise_i) == nothing
							break
						end
						stop -= 1
					end

					ylb[start:stop] .+= raise
					yub[start:stop] .+= raise
					first_decendant_drawn[parent] = true
				end
				y_bounds[sc,:,1] = ylb
				y_bounds[sc,:,2] = yub
				y_bounds[parent,:,1] = yub
			end
			area_plot(time_grid, yub, ybottom=ylb, color=GROUP_PALETTE[sc])
		end
	end
end


end # End of module
