
__precompile__()

module CopyNumberDeconvolution

export cn_deconvolution, cn_deconvolution_plot


# THIS CODE IS WORK-IN-PROGRESS 


# Math for expected lr, hsaf & maf values
expected_logratio(cf, copies, ploidy) = log2((cf * copies + ploidy - cf * ploidy) / ploidy)
expected_hsaf(cf, alt_cn, tot_cn) = (cf .* alt_cn .+ (1 .- cf)) ./ (cf .* tot_cn .+ (1 .- cf) .* 2)
balanced_hsaf_correction(depth) = 0.525 + (1.125 / depth) # simple empirically devised correction for expected hsaf level 0.5
function near_balanced_hsaf_correction(lvl, depth) # empirically devised correction for hsaf's close to 0.5
    balanced_correction = balanced_hsaf_correction(depth)
    to_correct_thresh = balanced_correction+depth/4000
    if lvl < to_correct_thresh
        return lvl + (to_correct_thresh - lvl) / 
                        (to_correct_thresh - 0.5) * 
                            (balanced_correction - 0.5)
    else
        return lvl
    end
end



# DECONVOLUTION

# Rules for model fit
# 1. Minimize abs(expected coverage - observed coverage)
# 2. Minimize abs(expected hsaf - observed hsaf)
# 3. Minimize number of CNAs needed for CN statuses
# 4. No large areas with 0 copies allowed
# 5. No allele specific copy number can get gain from 0

using IterTools: product
function subclone_and_allele_specific_cn_combinations(n_clones, max_cn)
    possible_tot_cns = collect(product(fill(0:max_cn, n_clones)...))[:];
    possible_cns = vcat(map(possible_tot_cns) do cns
        alts = map(c -> collect.(collect(product(c, 0:c))), cns)
        alt_combs = collect(product(map(c -> 1:length(0:c), cns)...))[:]
        permutedims(cat(dims=3, map(cmb -> hcat(map(i -> alts[i][cmb[i]], 1:length(alts))...), alt_combs)...), (3,2,1))
    end...);
end

function cn_deconvolution(cf::Vector{Float64}, dl::AbstractVector, depth::AbstractVector, 
                            segment_bounds, segment_vals,
                            sc_fracs::Array{Float64}, sc_parents::Vector{Int};
                            lr_const=1, hsaf_const=0.15, cna_const=0.05,
                            max_cn=7, max_0_seg_length=2e7)

    n_clones = size(sc_fracs, 1);
    n_segs = length(segment_bounds);

    # search space
    possible_cns = subclone_and_allele_specific_cn_combinations(n_clones, max_cn)

    # initialize result arrays
    cns = Array{Union{Nothing, Int}}(nothing, (n_segs, n_clones, 2))
    errs = fill(Inf, n_segs);

    # find best cn combination for each segment separately
    for seg_i in 1:n_segs

        # observed segment data
        len = segment_bounds[seg_i].stop - segment_bounds[seg_i].start;
        logr_obs = segment_vals[seg_i, :, 1];
        hsaf_obs = segment_vals[seg_i, :, 2];
        hsaf_obs_sd = sqrt.(segment_vals[seg_i, :, 4]);

        # keep track of best search result
        best_err = Inf;
        best_cns = nothing;

        # exhaustive search
        for i in 1:size(possible_cns, 1)

            tots_try = possible_cns[i,:,1]
            alts_try = possible_cns[i,:,2]

            if any(tots_try.==0) && len > max_0_seg_length; continue; end;

            # exp-obs lr penalty
            tots_exp = sum(sc_fracs .* tots_try, dims=1)[:];
            logr_exp = expected_logratio.(cf, tots_exp, 2) .+ dl;
            err = lr_const * maximum(abs.(logr_exp .- logr_obs));
            if err > best_err; continue; end;

            # exp-obs hsaf penalty
            alts_exp = sum(sc_fracs .* alts_try, dims=1)[:];
            hsaf_exp = expected_hsaf(cf, alts_exp, tots_exp) #(cf .* alts_exp .+ (1 .- cf).*1) ./ (cf .* tots_exp .+ (1 .- cf).*2);
            hsaf_exp = 0.5 .+ abs.(0.5 .- hsaf_exp);
            hsaf_exp = near_balanced_hsaf_correction.(hsaf_exp, depth)
            err += hsaf_const * maximum(abs.(hsaf_obs .- hsaf_exp) / hsaf_obs_sd);
            if err > best_err; continue; end;

            # minimum number of CNAs penalty
            min_cnas = 0
            for c in 1:n_clones
                if sc_parents[c]==0
                    continue # do not penalize pre-MRCA CNAs
                    # old_cns = (1,1)
                else
                    # discard if total cn is gained from 0
                    if tots_try[sc_parents[c]]==0 && tots_try[c]!=0; err=Inf; break; end

                    # discard if one of the two alleles is gained from 0
                    if (tots_try[sc_parents[c]]==alts_try[sc_parents[c]] || 
                        alts_try[sc_parents[c]]==0) && tots_try[c]!=alts_try[c]; err=Inf; break; end

                    # store parent's allele-specific copy numbers
                    old_cns = [alts_try[sc_parents[c]], tots_try[sc_parents[c]]-alts_try[sc_parents[c]]]
                end

                # store subclone's allele-specific copy numbers
                new_cns = [alts_try[c], tots_try[c]-alts_try[c]]

                # calculate minimum number of CNAs
                min_cnas += sum(abs.(new_cns .- old_cns))
            end
            err += cna_const * min_cnas;
            if err > best_err; continue; end;

            best_err = err;
            best_cns = possible_cns[i,:,:];
        end
        errs[seg_i] = best_err;
        cns[seg_i, :, :] .= best_cns;
    end

    return Array{Int}(cns), errs
end



# VISUALIZATION
d = readtsv("/home/joonatan/homo_sapiens/hg38.chrom.sizes")
const CHR_NAMES = d[:,1]
const CHR_SIZES = Int64.(d[:,2])
const CHR_STARTS = cumsum(CHR_SIZES) .- CHR_SIZES
const LEVEL_COLORS = Dict(0=>RGB(0, 0, 255), 
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
                    14=>RGB(0), 
                    999=>RGB(0, 255, 0),
                    998=>RGB(200, 0, 200),
                    997=>RGB(255, 255, 0),
                    996=>RGB(255, 153, 51),
                    995=>RGB(0),
                    994=>RGB(150)
)
const NOHIT_CN = 999
const BAD_DECONVOLUTION = 998

# todo to Plot2.jl
import PyCall
PyPlot = PyCall.pyimport("matplotlib.pyplot")
# function arrow_plot()
#     ...
# end

function cn_deconvolution_plot(out::String, samples::Vector{String},
                                cf::Vector{Float64}, dl::AbstractVector, depth::Vector{Int},
                                lr::Array{Float64}, lr_chr::Vector{String}, lr_pos::Vector{Int},
                                hsaf::Vector{Vector{Float64}}, hsaf_chr::Array{String}, hsaf_pos::Array{Int},
                                segment_bounds, cns;
                                max_cn=9, bad_deconvolutions=falses(size(cns, 1)),
                                plot_sex_chrs=false)

    n_samples = length(samples)
    n_segs = length(segment_bounds);
    n_clones = size(sc_fracs, 1);
    last_chr_i = plot_sex_chrs ? 24 : 22
    
    figure(out, size=(10, 2.7*(n_clones+n_samples))) do
        
        # Subclone allele-specific copy number profiles
        for sc_i in 1:n_clones
            # Raw copy numbers
            subplot(sc_i*2-1, (n_clones + n_samples)*2)  
                for (seg_i, s) in enumerate(segment_bounds)
                    if isnothing(cns[seg_i, 1, 1]); continue; end

                    # # dont draw segments with no differences between subclones
                    # if all(cns[seg_i, :, :] .== cns[seg_i, 1, :]'); continue; end;

                    cn = cns[seg_i, sc_i, 1]
                    color = LEVEL_COLORS[cn]

                    # draw segment
                    line_plot([s.chr_start+s.start, s.chr_start+s.stop], fill(cn, 2),
                                color=color, line_width=6)

                    # draw arrow indicating CNA from parent clone
                    if sc_parents[sc_i] != 0
                        middle = segment_bounds[seg_i].chr_start + 
                                (segment_bounds[seg_i].stop + segment_bounds[seg_i].start) / 2
                        cn_parent = cns[seg_i, sc_parents[sc_i], 1]
                        cn_diff = cns[seg_i, sc_i, 1] - cn_parent
                        cn_diff = cn_diff - sign(cn_diff)*0.4
                        PyPlot.arrow(middle, cn_parent, 0, cn_diff, width=4e6, 
                                        head_width=1.2e7, head_length=0.3, 
                                        color=hex(LEVEL_COLORS[cn_parent]), 
                                        zorder=10, length_includes_head=true)
                    end

                    # mark bad deconvolutions with "x"
                    if bad_deconvolutions[seg_i]
                        middle = ((s.chr_start+s.start)+(s.chr_start+s.stop))/2
                        scatter_plot([middle], [cn], size=15, marker="x", 
                                        color=LEVEL_COLORS[BAD_DECONVOLUTION], zorder=10)
                    end
                end

                genome_plot_config(CHR_NAMES[1:last_chr_i], CHR_SIZES[1:last_chr_i])
                grid(false, axis="y")
                remove_xtick_labels()
                yticks([0,1,2,3,4,5,6]) # TODO: change depending on detected levels
                ylim(-1,max_cn)
                xlim(0, CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i])
                ylabel("Copynumber")
                title("Clone "*('A':'Z')[sc_in])

            # HSAF levels
            subplot(sc_i*2, (n_clones + n_samples)*2)
                # store observed levels for later labeling
                observed_levels = Vector{NamedTuple{(:level, :cn_tot, :cn_alt), Float64, Int, Int}}()
                for (seg_i, s) in enumerate(segment_bounds)
                    if isnothing(cns[seg_i, 1, 1]); continue; end

                    # # dont draw segments with no differences between subclones
                    # if all(cns[seg_i, :, :] .== cns[seg_i, 1, :]'); continue; end;

                    cn_tot = cns[seg_i, sc_i, 1]
                    cn_alt = cns[seg_i, sc_i, 2]
                    color=LEVEL_COLORS[cn_tot]
                    lvl = expected_hsaf(1.0, cn_alt, cn_tot)
                    push!(observed_levels, (level=lvl, cn_tot=cn_tot, cn_alt=cn_alt))

                    # draw segment
                    line_plot([s.chr_start+s.start, s.chr_start+s.stop], fill(lvl, 2),
                                color=color, line_width=6)

                    # draw arrow indicating CNA from parent clone
                    if sc_parents[sc_i]!=0
                        middle = segment_bounds[seg_i].chr_start + 
                                    (segment_bounds[seg_i].stop + segment_bounds[seg_i].start) / 2
                        parent_tot = cns[seg_i, sc_parents[sc_i], 1]
                        parent_alt = cns[seg_i, sc_parents[sc_i], 2]
                        parent_lvl = expected_hsaf(1.0, parent_alt, parent_tot)
                        # push!(observed_levels, (level=parent_lvl, cn_tot=parent_tot, 
                        #                            cn_alt=parent_alt))
                        lv_diff = lvl - parent_lvl
                        lv_diff = lv_diff - sign(lv_diff)*0.03
                        PyPlot.arrow(middle, parent_lvl, 0, lv_diff, width=4e6, 
                                        head_width=1.2e7, head_length=0.02, 
                                        color=hex(LEVEL_COLORS[parent_tot]), 
                                        zorder=10, length_includes_head=true)
                    end

                    # mark bad deconvolutions with "x"
                    if bad_deconvolutions[seg_i]
                        middle = ((s.chr_start+s.start)+(s.chr_start+s.stop))/2
                        scatter_plot([middle], [lvl], marker="x", size=15, 
                                color=LEVEL_COLORS[BAD_DECONVOLUTION], zorder=10)
                    end
                end

                # right side y-axis labels for observed copy number states
                for l in unique(observed_levels)
                    label = "$(l.cn_alt)/$(l.cn_tot)"
                    text(CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i]+3e7, l.level, 
                            label, h_align="left", size=5, color=LEVEL_COLORS[l.cn_tot])
                end

                genome_plot_config(CHR_NAMES[1:last_chr_i], CHR_SIZES[1:last_chr_i])
                remove_xtick_labels()
                ylabel("HSAF")
                xlim(0, CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i])
                grid(false, axis="y")
        end

        # Sample LR scatters and expected segment levels
        for (si, sample) in enumerate(samples)
            # Coverage
            subplot(2*n_clones+2*si-1, (n_clones + n_samples)*2)
                # coverage logratio scatterplot
                genome_scatter_plot(lr_chr, lr_pos, lr[:, si])

                for (seg_i, s) in enumerate(segment_bounds)
                    if isnothing(cns[seg_i, 1, 1]); continue; end
                    
                    # if there are no CNAs between subclones, set segment color according to cn
                    if all(cns[seg_i, :, :] .== cns[seg_i, 1, :]')
                        cn = cns[seg_i, 1, 1]
                        color=LEVEL_COLORS[cn]
                    # otherwise set color to green
                    else
                        cn = sum(sc_fracs[:,si] .* cns[seg_i, :, 1])
                        color = LEVEL_COLORS[NOHIT_CN]
                    end

                    # draw expected segment level
                    line_plot([s.chr_start+s.start, s.chr_start+s.stop], 
                                fill(expected_logratio(cf[si], cn)+dl[si], 2), 
                                color=color, line_width=2)

                    # mark bad deconvolutions with "x"
                    if bad_deconvolutions[seg_i]
                        middle = ((s.chr_start+s.start)+(s.chr_start+s.stop))/2
                        scatter_plot([middle], [cn], size=15, marker="x", zorder=10, 
                                        color=LEVEL_COLORS[BAD_DECONVOLUTION])
                    end
                end

                # define ymin & ymax
                lvl_below_min_obs_lvl = expected_logratio(cf[si], max(0, minimum(cns[:,:,1])-1))+dl[si]
                lvl_above_max_obs_lvl = expected_logratio(cf[si], max(0, maximum(cns[:,:,1])+1))+dl[si]
                sex_chr_filter = plot_sex_chrs ? trues(length(lr_chr)) : !occursin.(r"^.+(X|Y)", lr_chr)
                min_lr_point = minimum(lr[sex_chr_filter, si])
                max_lr_point = maximum(lr[sex_chr_filter, si])
                ymin = max(lvl_below_min_obs_lvl, min_lr_point)
                ymax = min(lvl_above_max_obs_lvl, max_lr_point)
                
                # right side y-axis labels for expected copy number levels
                for l in cn_lvls[si]
                    if !(ymin < l.level < ymax); continue; end
                    # line_plot([0, CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i]], fill(l.level, 2), 
                    #            color=LEVEL_COLORS[l.cn], line_width=0.5)
                    text(CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i]]+3e7, l.level, l.cn, 
                            h_align="left", size=7, color=LEVEL_COLORS[l.cn])
                end

                ylim(ymin, ymax)
                xlim(0, CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i])
                ylabel("Coverage logratio")
                t = sample * join([", $(Int(round(sc_fracs[p_i, si]*100)))% $('@'+p_i)" 
                                    for p_i in 1:n_clones])
                title(t)

            # HSAF
            subplot(2*n_clones+2*si, (n_clones + n_samples)*2)
                # hsaf scatter plot
                genome_scatter_plot(hsaf_chr[si], hsaf_pos[si], hsaf[si])

                for (seg_i, s) in enumerate(segment_bounds)
                    if isnothing(cns[seg_i, 1, 1]); continue; end

                    # find expected levels
                    if all(cns[seg_i, :, :] .== cns[seg_i, 1, :]') # if there are no CNAs between subclones
                        cn_tot = cns[seg_i, 1, 1]
                        cn_alt = cns[seg_i, 1, 2]
                        lvl = expected_hsaf(cf[si], cn_alt, cn_tot)
                        # set color according to cn
                        color=LEVEL_COLORS[cn_tot]
                    else
                        cn_tot = cns[seg_i, :, 1]
                        cn_alt = cns[seg_i, :, 2]
                        cn_tot_avg = sum(sc_fracs[:,si] .* cn_tot);
                        cn_alt_avg = sum(sc_fracs[:,si] .* cn_alt);
                        lvl = expected_hsaf(cf[si], cn_alt_avg, cn_tot_avg);
                        lvl = 0.5 + abs(0.5 - lvl);
                        lvl = hsaf_level_correction(lvl, depth[si])
                        # set color to green
                        color = LEVEL_COLORS[NOHIT_CN]
                    end

                    # draw expected segment level
                    line_plot([s.chr_start+s.start, s.chr_start+s.stop], fill(lvl, 2), 
                                color=color, line_width=2)
                    
                    # mark bad deconvolutions with "x"
                    if bad_deconvolutions[seg_i]
                        middle = ((s.chr_start+s.start)+(s.chr_start+s.stop))/2
                        scatter_plot([middle], [lvl], size=15, marker="x", 
                                        color=LEVEL_COLORS[BAD_DECONVOLUTION], zorder=10)
                    end
                end

                ylabel("HSAF")
                xlim(0, CHR_STARTS[last_chr_i]+CHR_SIZES[last_chr_i])
        end
    end
end


end  # End of module
