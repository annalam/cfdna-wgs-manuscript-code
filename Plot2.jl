
__precompile__()

module Plot2

export figure, RGB, lighter_color
export line_plot, scatter_plot, bar_plot, violin_plot
export area_plot, matrix_plot
export beeswarm_plot, survival_arrow_plot, stacked_bar_plot, box_plot
export rectangle, arc, ellipse, errorbar, text
export xlim, ylim, xticks, yticks, subplot, ylabel, xlabel, title, grid
export genome_scatter_plot, genome_segments_plot, genome_plot_config
export suptitle, kde_plot,  genome_kde_plot, genome_plot_annotation, get_ylim, get_xlim

using Helpers, Printf, Survival, DelimitedFiles, CopyNum
import PyCall

global PyPlot = nothing;
global drawing_figure = false;

function figure(block::Function, path::String; size=(5, 3))
	global drawing_figure, PyPlot;
	if drawing_figure == true
		block()
	else
		path = expanduser(path)
		if !isdir(dirname(abspath(path)))
			error("Cannot create figure $(path). Directory does not exist.")
		end

		# Import matplotlib.pyplot if we have not done so yet
		if PyPlot == nothing
			#matplotlib = PyCall.pyimport("matplotlib")
			PyPlot = PyCall.pyimport("matplotlib.pyplot")
		end

		drawing_figure = true

		pyplot_fig = PyPlot.figure(figsize=size)
		try
			block()
			PyPlot.tight_layout()
			PyPlot.savefig(path)
			println(stderr, "Figure saved at $(path)")
		finally
			PyPlot.close()
			drawing_figure = false
		end
	end
	return nothing
end

figure(block::Function; kwargs...) = figure(block, expanduser("~/plot.pdf"); kwargs...)

xlabel(label::String) = PyPlot.xlabel(label)
ylabel(label::String) = PyPlot.ylabel(label)
function title(text::String; offset=0, location="center", size="large")
	PyPlot.title(text, pad=offset, loc=location, size=size)
end

function suptitle(text::String; fontsize=15)
	PyPlot.suptitle(text, fontsize=fontsize)
end

function xticks(ticks::AbstractVector; rotation=90, spine=true, grid=false, labels=nothing)
	if all(x -> x isa AbstractString, ticks)
		PyPlot.xticks(1:length(ticks), ticks, rotation=rotation)
		PyPlot.tick_params(axis="x", length=0)
	else
		PyPlot.xticks(ticks)
	end
	PyPlot.gca().spines["bottom"].set_visible(spine)
	PyPlot.gca().spines["top"].set_visible(false)
	if grid
		PyPlot.grid(true, axis="x", linestyle="dashed")
	else
		PyPlot.grid(false, axis="x")
	end
	if labels != nothing
		PyPlot.gca().set_xticklabels(labels)
	end
end

function yticks(ticks::AbstractVector; spine=true, grid=false, labels=nothing)
	if all(x -> x isa AbstractString, ticks)
		PyPlot.yticks(1:length(ticks), ticks)
		PyPlot.tick_params(axis="y", length=0)
	else
		PyPlot.yticks(ticks)
	end
	PyPlot.gca().spines["left"].set_visible(spine)
	PyPlot.gca().spines["right"].set_visible(false)
	if grid
		PyPlot.grid(true, axis="y", linestyle="dashed")
	else
		PyPlot.grid(false, axis="y")
	end
	if labels != nothing
		PyPlot.gca().set_yticklabels(labels)
	end
end

#remove_xtick_labels() = PyPlot.gca().set_xticklabels([]); PyPlot.gca().set_xticklabels([], minor=true)
#remove_ytick_labels() = PyPlot.gca().set_yticklabels([]); PyPlot.gca().set_yticklabels([], minor=true)

grid(args...; kwargs...) = PyPlot.grid(args...; kwargs...)

function xlim(low::Real, high::Real; log=NaN)
	if isfinite(low); PyPlot.xlim(left=low); end
	if isfinite(high); PyPlot.xlim(right=high); end
	if log > 1; PyPlot.semilogx(base=log); end
end

function ylim(low::Real, high::Real; log=NaN)
	if isfinite(low); PyPlot.ylim(bottom=low); end
	if isfinite(high); PyPlot.ylim(top=high); end
	if log > 1; PyPlot.semilogy(base=log); end
end
get_ylim() = PyPlot.gca().get_ylim()
get_ylim(axis) = PyPlot.gca().get_ylim()[axis]
get_xlim() = PyPlot.gca().get_xlim()
get_xlim(axis) = PyPlot.gca().get_xlim()[axis]

function subplot(panel::Integer, rows::Integer)
	PyPlot.subplot(rows, 1, panel)
end

function subplot(panel::Integer; rows::Integer=1, cols::Integer=1)
	PyPlot.subplot(rows, cols, panel)
end

struct RGB; r::UInt8; g::UInt8; b::UInt8; end
RGB(gray::Integer) = RGB(gray, gray, gray)
Base.broadcastable(r::RGB) = Ref(r)
hex(c::RGB) = @sprintf("#%02X%02X%02X", c.r, c.g, c.b)
hex(colors::Vector{RGB}) = hex.(colors)
function lighter_color(rgb::RGB; factor=0.2)
	RGB(min( round(rgb.r+(255-rgb.r)*factor, digits=0), 255 ), 
		min( round(rgb.g+(255-rgb.g)*factor, digits=0), 255 ), 
		min( round(rgb.b+(255-rgb.b)*factor, digits=0), 255 ))
end

function text(x, y, str; size=10, h_align="center", v_align="center", color=RGB(0), rotation="horizontal")
	PyPlot.text(x, y, str, fontsize=size, horizontalalignment=h_align, verticalalignment=v_align, color=hex(color), rotation=rotation)
end

function line_plot(x::AbstractVector, y::AbstractVector; color=RGB(0,0,0), line_width=1, capstyle="butt", style="solid")
	figure() do
		PyPlot.plot(x, y, color=hex(color), linewidth=line_width, solid_capstyle=capstyle, linestyle=style)
		PyPlot.grid(true, axis="both", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)
	end
end

function errorbar(args...; kwargs...)
	figure() do
		PyPlot.errorbar(args...; kwargs...)
	end
end

function genome_plot_config(chr_names, chr_lens)
	PyPlot.grid(true, axis="both", which="major", linestyle="dashed")
	PyPlot.xlim(0, sum(chr_lens))
	PyPlot.tick_params(axis="x", which="both", length=0)
	PyPlot.gca().set_xticks(cumsum(chr_lens))
	PyPlot.gca().set_xticklabels([])
	PyPlot.gca().set_xticks(cumsum(chr_lens) .- chr_lens ./ 2, minor=true)
	PyPlot.gca().set_xticklabels(
		[r"(19|21|M)$" in chr ? "" : replace(chr, "chr", "")
		for chr in chr_names], minor=true)
end

function genome_scatter_plot(chromosome::Vector{String},
	position::AbstractVector, value::AbstractVector;
	chr_sizes="~/homo_sapiens/hg38.chrom.sizes", color=RGB(0), line=false)

	@assert(length(chromosome) == length(position))
	d = readdlm(expanduser(chr_sizes))
	chr_names = d[:, 1]
	chr_lens = Int32.(d[:, 2])
	chr_start = Dict(chr_names[c] => sum(chr_lens[1:c-1])
		for c in 1:length(chr_names))

	gpos = [chr_start[chromosome[k]] + position[k]
			for k in 1:length(chromosome)]

	figure() do
		if line
			order = sortperm(gpos)
			line_plot(gpos[order], value[order], color=color)
		else
			scatter_plot(gpos, value, color=color)
		end
		genome_plot_config(chr_names, chr_lens)
	end
end

function genome_segments_plot(segments::Vector{Segment},
	seg_values::Vector; chr_sizes="~/homo_sapiens/hg38.chrom.sizes",
	color=RGB(255, 0, 0), line_width=2)

	@assert(length(segments) == length(seg_values))
	d = readdlm(expanduser(chr_sizes))
	chr_names = d[:, 1]
	chr_lens = Int32.(d[:, 2])
	chr_start = Dict(chr_names[c] => sum(chr_lens[1:c-1])
		for c in 1:length(chr_names))

	figure() do
		for (seg, value) in zip(segments, seg_values)
			start = chr_start[seg.chromosome] + seg.start
			stop = chr_start[seg.chromosome] + seg.stop
			line_plot([start, stop], fill(value, 2),
				color=color, line_width=line_width)
		end
		genome_plot_config(chr_names, chr_lens)
	end
end

using Distributions
function genome_kde_plot(chromosome::AbstractVector, position::AbstractVector;
	chr_sizes="~/homo_sapiens/hg38.chrom.sizes", sigma=1e6)

	@assert(length(chromosome) == length(position))
	d = readdlm(expanduser(chr_sizes))
	chr_names = d[:, 1]
	chr_lens = Int32.(d[:, 2])
	chr_start = map(c -> sum(chr_lens[1:c-1]), eachindex(chr_lens))

	gpos = map(1:length(chromosome)) do k
		chr = findone(chr_names, chromosome[k])
		start = chr_start[chr]
		pos = start + position[k]
		stop = start + chr_lens[chr]
		return [pos, start, stop]
	end
	
	x = 0:1e6:(chr_start[end]+chr_lens[end])
	y = zeros(length(x))
	for gp in gpos
		i = gp[2] .< x .< gp[3]
		yy = pdf(Normal(gp[1], sigma), x[i])
		y[i] += yy
	end

	figure() do
		line_plot(x, y)
		genome_plot_config(chr_names, chr_lens)
	end
end

function genome_plot_annotation(chromosome::AbstractVector, position::AbstractVector, annotation::AbstractVector; 
	chr_sizes="~/homo_sapiens/hg38.chrom.sizes", color=RGB(255,0,0), fontsize=10, line_width=1, rotation="horizontal")

	@assert(length(chromosome) == length(position))
	d = readdlm(expanduser(chr_sizes))
	chr_names = d[:, 1]
	chr_lens = Int32.(d[:, 2])
	chr_start = Dict(chr_names[c] => sum(chr_lens[1:c-1])
		for c in 1:length(chr_names))

	gpos = [chr_start[chromosome[k]] + position[k]
			for k in 1:length(chromosome)]
	
	ylim = [l for l in PyPlot.gca().get_ylim()]

	for (i, pos) in enumerate(gpos)
		line_plot([pos, pos], ylim, color=color, line_width=line_width)
		text(pos, ylim[2], annotation[i], size=fontsize, v_align="bottom", rotation=rotation)
	end
end

using KernelDensity
function kde_plot(x::AbstractArray; bandwidth=0.5, line_width=1, npoints=2048, color=RGB(0))
	U = kde(Float64.(x), bandwidth=bandwidth, npoints=npoints)
	line_plot(U.x, U.density, line_width=line_width, color=color)
end

function scatter_plot(x::AbstractVector, y::AbstractVector;
	size=2, linewidths=1.0, color=RGB(0,0,0), opacity=1.0, border_color=nothing, fill_color=nothing, marker="o", 
	xerror=nothing, yerror=nothing, error_line_width=1, error_color=RGB(0,0,0), grid=true, zorder=1)

	@assert(length(x) == length(y))

	# This is a workaround for the fact that Matplotlib misindexes the
	# "error_color" vector if some error ranges are NaN.
	if xerror != nothing
		@assert(xerror isa AbstractMatrix && Base.size(xerror, 1) == length(x))
		for k in 1:Base.size(xerror, 1)
			if any(isnan, xerror[k, :]); xerror[k, :] .= 0; end
		end
	end
	if yerror != nothing
		@assert(yerror isa AbstractMatrix && Base.size(yerror, 1) == length(x))
		for k in 1:Base.size(yerror, 1)
			if any(isnan, yerror[k, :]); yerror[k, :] .= 0; end
		end
	end

	figure() do
		if xerror isa AbstractMatrix
			PyPlot.hlines(y, xmin=xerror[:, 1], xmax=xerror[:, 2], colors=hex(error_color), linewidth=error_line_width, zorder=0)
		end
		if yerror isa AbstractMatrix
			PyPlot.vlines(x, ymin=yerror[:, 1], ymax=yerror[:, 2], colors=hex(error_color), linewidth=error_line_width, zorder=0)
		end

		if grid == true
			PyPlot.grid(true, axis="both", linestyle="dashed")
		else
			PyPlot.grid(false)
		end
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		edgecolors = all(isa.(border_color, RGB)) ? hex(border_color) : "none"
		facecolors = isnothing(fill_color) ? hex(color) : isa(fill_color, RGB) ? hex(fill_color) : fill_color
		color = facecolors=="none" ? "none" : hex(color)
		PyPlot.scatter(x, y, s=size, linewidths=linewidths, marker=marker, c=color,
			alpha=opacity, edgecolors=edgecolors, facecolors=facecolors, zorder=zorder)
	end
end

function bar_plot(values::AbstractVector; box_width=0.8, color=RGB(0,0,0))
	figure() do
		PyPlot.bar(1:length(values), values, color=hex(color),
			width=box_width, linewidth=0,
			bottom=0.0000000000000001, zorder=10)
		PyPlot.xlim(0, length(values) + 1)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)
	end
end

function stacked_bar_plot(values::AbstractMatrix; box_width=0.8, colors=[])
	B = size(values, 1)    # Number of bars
	L = size(values, 2)    # Number of levels

	if isempty(colors)
		colors = [RGB(round(Int, (l - 1) / L) * 255) for l in 1:L]
	end

	@assert(colors isa AbstractVector{RGB})
	if length(colors) != L
		error("Bar plot has $L levels but $(length(colors)) colors were provided in 'colors' keyword argument.")
	end

	figure() do
		bottoms = zeros(B)
		PyPlot.xlim(0, B + 1)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		for level in 1:L
			PyPlot.bar(1:B, values[:, level], bottom=bottoms, color=hex(colors[level]), zorder=10)
			bottoms .+= values[:, level]
		end
	end
end

function box_plot(groups...; show_means=false, show_outliers=true)
	figure() do
		PyPlot.boxplot(groups, showmeans=show_means, showfliers=show_outliers)
	end
end

function area_plot(x::AbstractVector, y::AbstractVector; xbottom=nothing, ybottom=nothing, color=RGB(0))

	figure() do
		if xbottom == nothing && ybottom == nothing
			PyPlot.fill_between(x, y, facecolor=hex(color))
		end
		if xbottom isa Number || xbottom isa AbstractVector
			PyPlot.fill_betweenx(y, xbottom, x, facecolor=hex(color))
		end
		if ybottom isa Number || ybottom isa AbstractVector
			PyPlot.fill_between(x, ybottom, y, facecolor=hex(color))
		end
	end
end

#function histogram_plot(values::Array; path="~/plot.pdf")
	#centers =
	#_, h = hist(values, edges); centers = edges[1:end-1] + step(edges) / 2
	#histogram_plot(values, path=path, box_width=1.0)
#end

function violin_plot(bins::AbstractVector, density::AbstractArray; colors=[], halved=false)

	V = size(density, 2);   # How many violins?

	bins = Float64.(bins); density = Float64.(density)
	for v in 1:V
		density[:, v] /= maximum(density[:, v]) * 2.05
	end

	if isempty(colors)
		colors = map(v -> RGB(0, 0, 0), 1:V)
	end

	figure() do
		PyPlot.xlim(0.25, V + 0.75)
		PyPlot.tick_params(axis="x", length=0)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		for v in 1:V
			if halved
				x = vcat(v .- density[:, v], v, v, v .- density[1, v])
				y = vcat(bins, bins[end], bins[1], bins[1])
			else
				x = vcat(v .- density[:, v], v .+ density[end:-1:1, v])
				y = vcat(bins, bins[end:-1:1])
			end
			PyPlot.fill(x, y, hex(colors[v]), zorder=10)
		end
	end
end

function rectangle(x1::Real, x2::Real, y1::Real, y2::Real; edge_width=0, color=RGB(0, 0, 0), opacity=1.0)
	width = abs(x2 - x1); height = abs(y2 - y1);
	x = min(x1, x2); y = min(y1, y2);
	rect = PyPlot.matplotlib.patches.Rectangle((x, y), width, height, linewidth=edge_width, facecolor=hex(color), alpha=opacity)
	PyPlot.gca().add_patch(rect)
end

function arc(x1::Real, x2::Real, y1::Real, y2::Real; color=RGB(0))
	@assert(y2 >= y1)
	patch = PyPlot.matplotlib.patches.Arc(((x1 + x2) / 2, y1), abs(x2 - x1), y2 - y1, theta2=180)
	PyPlot.gca().add_patch(patch)
end

function ellipse(x::Real, y::Real, rx::Real, ry::Real; color=RGB(0, 0, 0))
	patch = PyPlot.matplotlib.patches.Ellipse([x, y], 2*rx, 2*ry,
		facecolor=hex(color))
	PyPlot.gca().add_patch(patch)
end

function beeswarm_plot(groups...; point_size=2, color=RGB(0,0,0))
	x = zeros(0); y = zeros(0);
	for (group, values) in enumerate(groups)
		x = vcat(x, group .+ clamp.(randn(length(values)) / 10, -0.4, 0.4))
		y = vcat(y, values)
	end
	figure() do
		PyPlot.tick_params(axis="x", length=0)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)
		PyPlot.scatter(x, y, s=point_size, c=hex(color), marker=".", zorder=10)
	end
end

function survival_arrow_plot(survival::Vector{SurvivalTime}; color=RGB(0,0,0), start_time=[], head_length=NaN)

	S = length(survival)
	max_time = maximum(abs.(survival))
	if isnan(head_length); head_length = max_time / 20; end

	if isa(color, RGB); color = fill(color, S); end

	if isempty(start_time)
		start_time = zeros(S)
	end

	figure() do
		PyPlot.xlim(0.5, S + 0.5)
		PyPlot.ylim(0, max_time * 1.1)
		PyPlot.tick_params(axis="x", length=0)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		for s in 1:S
			if survival[s] == event(0); continue; end
			PyPlot.arrow(s, start_time[s], 0, abs(survival[s]),
				color=hex(color[s]), width=0.3,
				head_width=(is_censored(survival[s]) ? 0.8 : 0),
				head_length=(is_censored(survival[s]) ? head_length : 0))
		end
	end
end


export MatrixLayer, matrix_plot
export GLYPH_NONE, GLYPH_TILE, GLYPH_SQUARE, GLYPH_SQUARE_BELOW, GLYPH_SQUARE_ABOVE, GLYPH_HORIZONTAL_BOX

const GLYPH_NONE = 0
const GLYPH_TILE = 1
const GLYPH_SQUARE = 2
const GLYPH_SQUARE_BELOW = 3
const GLYPH_SQUARE_ABOVE = 4
const GLYPH_HORIZONTAL_BOX = 5

struct MatrixLayer
	glyph::Array{Int8}
	color::Array{RGB}
end
MatrixLayer(rows::Integer, cols::Integer; glyph=GLYPH_NONE) =
	MatrixLayer(fill(glyph, rows, cols), fill(RGB(255, 255, 255), rows, cols))
Base.getindex(layer::MatrixLayer, rows, cols) = MatrixLayer(layer.glyph[rows, cols], layer.color[rows, cols])

function render_glyph(svg::IO, glyph::Int8, color::RGB,
	x::Int64, y::Int64, cell_width::Int64, cell_height::Int64)

	rgb = "rgb($(color.r), $(color.g), $(color.b))"
	if glyph == GLYPH_TILE
		write(svg, """<rect x="$(x)" y="$(y)" width="$(cell_width)" height="$(cell_height)" style="stroke: rgb(255, 255, 255); fill: $rgb" />""")
	elseif glyph == GLYPH_SQUARE
		write(svg, """<rect x="$(x + cell_width/4)" y="$(y + cell_height/4)" width="$(cell_width/2)" height="$(cell_height/2)" style="stroke: none; fill: $rgb" />""")
	elseif glyph == GLYPH_SQUARE_BELOW
		write(svg, """<rect x="$(x + cell_width/6)" y="$(y + cell_height/6)" width="$(cell_width/2)" height="$(cell_height/2)" style="stroke: none; fill: $rgb" />""")
	elseif glyph == GLYPH_SQUARE_ABOVE
		write(svg, """<rect x="$(x + cell_width*2/6)" y="$(y + cell_height*2/6))" width="$(cell_width/2)" height="$(cell_height/2)" style="stroke: none; fill: $rgb" />""")
	elseif glyph == GLYPH_HORIZONTAL_BOX
		write(svg, """<rect x="$(x)" y="$(y + cell_height/4)" width="$(cell_width)" height="$(cell_height/2)" style="stroke: none; fill: $rgb" />""")
	end
end

function matrix_plot(layers::Array{MatrixLayer}; path="~/plot.svg", cell_width=50, cell_height=50)
	R, C = size(layers[1].glyph)
	canvas_w = C * cell_width
	canvas_h = R * cell_height
	svg = open(expanduser(path), "w")
	write(svg, """<svg width="$(canvas_w)" height="$(canvas_h)">""")
	for layer in layers
		@assert(size(layer.glyph) == (R, C))
		for r in 1:R
			for c in 1:C
				x = (c - 1) * cell_width
				y = (r - 1) * cell_height
				render_glyph(svg, layer.glyph[r, c], layer.color[r, c], x, y, cell_width, cell_height)
			end
		end
	end
	write(svg, "</svg>")
	close(svg)
end

matrix_plot(layer::MatrixLayer; kwargs...) = matrix_plot([layer]; kwargs...)

end
