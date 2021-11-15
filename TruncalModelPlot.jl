
__precompile__()

module TruncalModelPlot

export truncal_model_plot


using Helpers, Printf, CopyNum, Plot2



# Math for expected lr, hsaf & maf values
function raw_cn(lr, cf, dl; sex_chr=false)
    if sex_chr
        (2 .^ (lr .- dl) .- 1 + cf) ./ cf
    else
        (2 .^ (lr .- dl .+ 1) .- 2 * (1-cf)) ./ cf
    end
end
expected_logratio(ctfr, copies, ploidy) = log2((ctfr * copies + ploidy - ctfr * ploidy) / ploidy)
expected_hsaf(cf, alt_cn, tot_cn) = (cf*alt_cn+(1-cf)*1) / ( cf*(tot_cn) + (1-cf)*2 )
expected_maf(cf, alt_cn, tot_cn) = (cf*alt_cn) / ( cf*(tot_cn) + (1-cf)*2 )
balanced_hsaf_correction(depth) = 0.525 + (1.125 / depth) # simple empirically devised correction for expected hsaf level 0.5



# Chromosome coordinates
d = readtsv("/home/joonatan/homo_sapiens/hg38.chrom.sizes")
const CHR_NAMES = d[:,1]
const CHR_SIZES = Int64.(d[:,2])
const CHR_STARTS = cumsum(CHR_SIZES) .- CHR_SIZES



# Copy number colors
const LEVEL_COLORS = Dict(
                    nothing=>RGB(0, 255, 0),
                    0=>RGB(0, 0, 255), 
                    1=>RGB(120, 120, 255), 
                    2=>RGB(100), 
                    3=>RGB(255, 180, 180), 
                    4=>RGB(255, 120, 120), 
                    5=>RGB(255, 0, 0), 
                    6=>RGB(180, 0, 0), 
                    7=>RGB(120, 0, 0), 
                    8=>RGB(60, 0, 0), 
                    9=>RGB(40, 0, 0), 
                    10=>RGB(20, 0, 0),
                    11=>RGB(10, 0, 0), 
                    12=>RGB(0), 
                    13=>RGB(0), 
                    14=>RGB(0)
)



# Assigning segment and mutations copy numbers

# todo: move to CopyNum
# Check what expected allele-specific copy number state values each segment matches, if any
using Distributions: Normal, quantile
function segment_copynumbers(segment_vals::Array{Float64, 3}, cn_lvls, hsaf_lvls; 
                                lr_q_thresh=0.2, hsaf_q_thresh=0.25, adj_lr_lvl_thresh=0.25)
    
    # segment values are given in a 3D array where dimensions are:
        # 1. segments
        # 2. samples
        # 3. lr mean, hsaf mean, lr var, hsaf var

    # initialize return arrays
    tot_cns = Array{Union{Nothing, Int}}(
        nothing, (size(segment_vals, 1), size(segment_vals, 2)))
    alt_cns = copy(tot_cns)

    for seg_i in 1:size(segment_vals, 1)
        for si in 1:size(segment_vals, 2) # si = sample index
            vals = segment_vals[seg_i, si, :]
            if any(isnan.(vals)); continue; end

            # segment lr values at quantile thresholds
            lr_quantiles = quantile(Normal(vals[1], sqrt(vals[3])), [lr_q_thresh, 1-lr_q_thresh])
            # expected lr levels inside threholds
            lr_cn_valid = map(enumerate(cn_lvls[si])) do (i, l)
                # the last level is just for checking if second last segments are valid
                if i == length(cn_lvls[si]); return false; end

                # check that the level is inside quantile thresholds
                distr_valid = lr_quantiles[1] .< l.level .< lr_quantiles[2]

                # check that adjacent levels are not too close
                upper_lvl_valid = vals[1] < l.level + (cn_lvls[si][i+1].level - l.level) * adj_lr_lvl_thresh
                if l.cn != 0
                    lower_lvl_valid = vals[1] > l.level - (l.level - cn_lvls[si][i-1].level) * adj_lr_lvl_thresh
                else
                    # if the level is 0, we can't check if the level below is too close
                    lower_lvl_valid = true
                end
                return distr_valid && upper_lvl_valid && lower_lvl_valid
            end
            lr_cn = cn_lvls[si][lr_cn_valid]

            # segment lr values at quantile thresholds
            saf_quantiles = quantile(Normal(vals[2], sqrt(vals[4])), [hsaf_q_thresh, 1-hsaf_q_thresh])
            # expected hsaf levels inside thresholds
            saf_cn = filter(l -> saf_quantiles[1] .< l.level .< saf_quantiles[2], hsaf_lvls[si]);

            # assign a final sigle allele-specific state for each segment
            tot_cn_hits = Int64[]
            alt_cn_hits = Int64[]
            for lr_hit in lr_cn, saf_hit in saf_cn
                if lr_hit.cn != saf_hit.cn; continue; end
                push!(tot_cn_hits, lr_hit.cn)
                push!(alt_cn_hits, saf_hit.alt_cn)
            end
            if length(tot_cn_hits) == 1
                tot_cns[seg_i, si] = tot_cn_hits[1]
                alt_cns[seg_i, si] = alt_cn_hits[1]
                continue
            else
                tot_cns[seg_i, si] = nothing
                alt_cns[seg_i, si] = nothing
                continue
            end
        end
    end

    return tot_cns, alt_cns
end
function segment_copynumbers(segment_vals::Array{Float64, 2}, cn_lvls, hsaf_lvls; 
                                lr_q_thresh=0.2, hsaf_q_thresh=0.25, adj_lr_lvl_thresh=0.25)
    # segment values are given in a 2D array where dimensions are:
        # 1. segments
        # 2. lr mean, hsaf mean, lr var, hsaf var
    
    # add "empty" dimensions
    if !(typeof(cn_lvls[1]) <: AbstractArray) && length(size(segment_vals)) == 2
        cn_lvls = [cn_lvls]
        hsaf_lvls = [hsaf_lvls]
        segment_vals = reshape(segment_vals, size(segment_vals, 1), 1, size(segment_vals, 2))
    end

    tot_cns, alt_cns = segment_copynumbers(segment_vals, cn_lvls, hsaf_lvls, 
                                            lr_q_thresh=lr_q_thresh, 
                                            hsaf_q_thresh=hsaf_q_thresh, 
                                            adj_lr_lvl_thresh=adj_lr_lvl_thresh)

    # remove the added dimension
    if size(segment_vals, 2) == 1
        tot_cns, alt_cns = tot_cns[:], alt_cns[:]
    end

    return tot_cns, alt_cns
end

# Find mutations in segments with each allele-specific copy number state
function muts_in_cns(chr, pos, segment_bounds, seg_tot_cns, seg_alt_cns)
    # segment_bounds needs to be sorted by genomic position
    indices = Int64[]
    tot_cns = zeros(Int64, 0, size(seg_tot_cns, 2))
    alt_cns = zeros(Int64, 0, size(seg_tot_cns, 2))

    for i in 1:length(chr), (seg_i, seg) in enumerate(segment_bounds)
        if CHR_STARTS[findone(CHR_NAMES, chr[i])] != seg.chr_start; continue; end
        if !(seg.start < pos[i] < seg.stop); continue; end
        tot_cn = seg_tot_cns[seg_i, :]
        alt_cn = seg_alt_cns[seg_i, :]
        if any(tot_cn .== nothing); continue; end
        push!(indices, i)
        tot_cns = vcat(tot_cns, transpose(tot_cn))
        alt_cns = vcat(alt_cns, transpose(alt_cn))
    end
    if size(tot_cns, 2) == 1
        tot_cns, alt_cns = tot_cns[:], alt_cns[:]
    end
    return indices, tot_cns, alt_cns
end



# Binomial kernel density estimate for visualizing MAFs
using Distributions: Binomial, pdf
using Interpolations: LinearInterpolation
function binomial_kde(grid, success_fractions::Vector{Float64}, trials::Vector{Int})
    density = zeros(length(grid))
    for i in 1:length(success_fractions)
        f = success_fractions[i]
        t = trials[i]
        if f==0 || t==0; continue; end
        p = pdf(Binomial(t, f), 0:t)
        interp = LinearInterpolation((0:t)./t, p)
        density += interp.(grid)
    end
    density = density ./ length(success_fractions) .* 100
    return density
end


function truncal_model_plot(out::String, cf::Float64, dl::Real, depth::Real,
                    lr::Vector{Float64}, lr_chr::Vector{String}, lr_pos::Vector{Int},
					hsaf::Vector{Float64}, hsaf_chr::Vector{String}, hsaf_pos::Vector{Int},
					maf::Vector{Float64}, m_tot::Vector{Int}, m_chr::Vector{String}, m_pos::Vector{Int}, # m_tot: number of mutant + reference reads
					segment_bounds, segment_vals; sample=nothing,
                    max_cn=9, adaptive_axis_limits=true, min_visualized_muts=100,
                    lr_q_thresh=0.2, hsaf_q_thresh=0.25, adj_lr_lvl_thresh=0.25,
					show_thresholds=false)

	# Calculate expected coverage logratio and HSAF levels, given cacer fraction and diploid level
    #max_cn+1 because assigning segment to copy number also looks at how far it is from adjacent levels
	cn_lvls = [(cn=copies, level=(expected_logratio(cf, copies, 2) + dl)) for copies in 0:max_cn+1];
    hsaf_lvls = NamedTuple{(:cn, :alt_cn, :level), Tuple{Int64, Int64, Float64}}[]
    for tot in 0:max_cn+1, alt in Int((tot + tot%2)/2):tot
        level = expected_hsaf(cf, alt, tot)
        if level==0.5
            level = balanced_hsaf_correction(depth)
        end
        push!(hsaf_lvls, (cn=tot, alt_cn=alt, level=level))
    end

	# Find which copy number state each segment represents
	s_cns, s_alt_cns = segment_copynumbers(segment_vals, cn_lvls, hsaf_lvls, 
        lr_q_thresh=lr_q_thresh, hsaf_q_thresh=hsaf_q_thresh, adj_lr_lvl_thresh=adj_lr_lvl_thresh);
	# All observed total copy numbers
	obs_cns	= unique(s_cns[s_cns.!=nothing]);

	# Find mutations from each total copy number
	m_idx, m_cns, _ = muts_in_cns(m_chr, m_pos, segment_bounds, s_cns, s_alt_cns);
	# Calculate expected MAF for each copy number state
    maf_lvls = NamedTuple{(:cn, :alt_cn, :level), Tuple{Int64, Int64, Float64}}[]
    for tot in unique(m_cns), alt in 1:tot
        push!(maf_lvls, (cn=tot, alt_cn=alt, level=expected_maf(cf, alt, tot)))
    end

	# Filter out large unobserved copy number states
	if length(obs_cns) != 0
		filter!(l -> l.cn <= max(4, maximum(obs_cns)), cn_lvls);
		filter!(l -> l.cn <= max(4, maximum(obs_cns)), hsaf_lvls);
		filter!(l -> l.cn in obs_cns, maf_lvls);
	end

	figure(out, size=(9,8)) do

		# Coverage logratio scatter plot
		subplot(1, 4)
            if adaptive_axis_limits
                ymin = floor(cn_lvls[2].level, digits=1) - floor(abs(dl - cn_lvls[2].level), digits=1) / 3
                ymax = floor(cn_lvls[end].level + 0.1, digits=1) + floor(abs(dl - cn_lvls[end].level) + 0.1, digits=1) / 3
            else
                ymin = minimum(lr)
                ymax = maximum(lr)
            end

			# Plot horizontal lines for expected levels
			for l in cn_lvls
				if l.level < ymin; continue; end
				line_plot([0, CHR_STARTS[25]], fill(l.level, 2), color=LEVEL_COLORS[l.cn], line_width=0.5)
				text(CHR_STARTS[25]+3e7, l.level, l.cn, h_align="left", size=7, color=LEVEL_COLORS[l.cn])
			end

			# Plot segment medians and thresholds
			for i in eachindex(segment_bounds)
				b = segment_bounds[i]
				if b.chr_start in CHR_STARTS[23:25]; continue; end;
				start, stop = b.chr_start .+ [b.start, b.stop]
				lr_mean, var = segment_vals[i, [1,3]]
				color=LEVEL_COLORS[s_cns[i]]
				line_plot([start, stop], fill(lr_mean, 2), color=color, line_width=2)
				if show_thresholds
					# Horizontal line showing the threshold quantiles of the segment
					q1, q2 = quantile(Normal(lr_mean, sqrt(var)), lr_q_thresh)
					line_plot(fill((start+stop)/2, 2), [q1, q2], color=color, line_width=1.2)
				end
			end

			# Scatter plot of coverage logratio values
			genome_scatter_plot(lr_chr, lr_pos, lr)

			grid(b=false, axis="y")
			ylim(ymin, ymax)
			ylabel("Coverage logratio")
            if sample==nothing
			    title(@sprintf("Cancer fraction: %d%%, Diploid level: %.2f", cf*100, dl))
            else
                title(@sprintf("%s - Cancer fraction: %d%%, Diploid level: %.2f", sample, cf*100, dl))
            end

		# HSAF scatter plot
		subplot(2, 4)
			ymax = round(maximum(hsaf) + 0.05, digits=1)

			# Plot horizontal lines for expected levels
			for l in hsaf_lvls
				if l.level > ymax; continue; end
				if l.cn > 4
                    if sum(s_cns .== l.cn) == 0
                        continue
                    elseif l.alt_cn > maximum(s_alt_cns[s_cns .== l.cn])
                        continue
                    end
                end
                allele_cns = [l.alt_cn, l.cn-l.alt_cn]
				label = "$(maximum(allele_cns))+$(minimum(allele_cns))"
				line_plot([0, CHR_STARTS[25]], fill(l.level, 2), color=LEVEL_COLORS[l.cn], line_width=0.5)
				text(CHR_STARTS[25]+3e7, l.level, label, h_align="left", size=5, color=LEVEL_COLORS[l.cn])
			end

			# Plot segment medians
			for i in eachindex(segment_bounds)
				b = segment_bounds[i]
				if b.chr_start in CHR_STARTS[23:25]; continue; end;
				start, stop = b.chr_start+b.start, b.chr_start+b.stop
				saf, var = segment_vals[i, [2,4]]
				if isnan(var); continue; end
				color=LEVEL_COLORS[s_cns[i]]
				line_plot([start, stop], fill(saf, 2), color=color, line_width=2)
				if show_thresholds
					# Horizontal line showing the threshold quantiles of the segment
					q1, q2 = quantile(Normal(saf, sqrt(var)), hsaf_q_thresh)
					line_plot(fill((start+stop)/2, 2), [q1, q2], color=color, line_width=1.2)
				end
			end

            # Scatter plot of hsaf values
			genome_scatter_plot(hsaf_chr, hsaf_pos, hsaf)

			grid(b=false, axis="y")
			ylim(0.5, ymax)
			ylabel("Hetz SNP allele fraction")
            xlabel("Chromosome")

		# HSAF + CN scatter plot
		subplot(3, rows=2, cols=2)
			# Plot expected levels as circles
            y = [l.level for l in hsaf_lvls]
            x = [cn_lvls[findone(lr_l -> lr_l.cn==l.cn, cn_lvls)].level for l in hsaf_lvls]
			x = raw_cn.(x, cf, dl) # convert coverage logratio to copy number scale
			c = [lighter_color(LEVEL_COLORS[l.cn]) for l in hsaf_lvls]
			scatter_plot(x, y, border_color=c, fill_color="none", size=150)

			# Plot segment points
			x = [raw_cn(lr, cf, dl) for lr in segment_vals[:, 1]] # raw segment CN values
			y = segment_vals[:, 2] # segment hsaf values
			s = [1+((s.stop-s.start) / 1e7) for s in segment_bounds] # point sizes scaled with segment lengths
			c = [LEVEL_COLORS[cn] for cn in s_cns]
			scatter_plot(x, y, color=c, size=s)
			grid(b=false)

			if length(obs_cns) != 0
				xlim(-0.2, maximum(obs_cns)+0.5)
			end
			ylim(0.5, maximum(segment_vals[:,2])+0.02)
			ylabel("Hetz SNP allele fraction")
			xlabel("Copynumber")

		# MAF distributions
		subplot(4, rows=2, cols=2)
			legend=""
            xmax, ymax = 0.1, 0
			for cn in obs_cns
				n_mut, n_seg = sum(m_cns.==cn), sum(s_cns.==cn)
				if n_mut < min_visualized_muts; continue; end;

				ii = m_idx[m_cns.==cn]
				f, tot = maf[ii], m_tot[ii]
                if maximum(f) > xmax
                    xmax = maximum(f)
                end
                x = 0:0.001:1
                y = binomial_kde(x, f, tot)
                y = y ./ maximum(y) .+ get_ylim()[2]
                line_plot(x, y, color=LEVEL_COLORS[cn], line_width=2)
                maf_lvls_cn = filter(l -> l.cn==cn, maf_lvls)
                for l in maf_lvls_cn
                    d = y[findfirst(x .> l.level)]
                    label = "$(l.alt_cn)/$(l.cn)"
                    if d-minimum(y) > 0.5
                        tick_end=d-0.2
                        v_align="top"
                    else
                        tick_end=d+0.2
                        v_align="bottom"
                    end
                    line_plot(fill(l.level ,2), [d, tick_end], color=LEVEL_COLORS[cn], line_width=2)
                    text(l.level, tick_end-0.05, label, size=6, v_align=v_align)
                end

                legend = "$n_seg segments\n$n_mut mutations"
                text(0, get_ylim()[2]-0.1, legend, size=6, 
                        v_align="top", h_align="left", color=LEVEL_COLORS[cn])
			end
			grid(b=false)
			xlim(0, xmax)
			ylim(0.9, NaN)
			yticks([])
            ylabel("Density");xlabel("MAF")
	end
end



end  # End of module
