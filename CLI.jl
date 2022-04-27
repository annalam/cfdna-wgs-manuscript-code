
__precompile__()

module CLI

using Helpers

export subcommands

struct Subcommand
	name::String
	#func::Function
	positional_names::Array{String}
	positional_types::Array{String}
	varargs_name::String
	varargs_type::String
	#options::Dict{String, String}  # FIXME: No introspection yet
end

function Subcommand(func::Function)
	methods = Base.methods(func).ms
	method = methods[1]
	@assert(length(methods) == 1 || PROGRAM_FILE != "")
	if length(methods) > 1
		matches = filter(m -> endswith(string(m.file), PROGRAM_FILE), methods)
		if length(matches) != 1
			error("Found multiple method matches:\n$(matches)")
		end
		method = matches[1]
	end

	slots = Base.arg_decl_parts(method)[2][2:end]
	positional_names = map(x -> x[1], slots)
	positional_types = map(x -> x[2], slots)
	varargs_name = ""; varargs_type = ""
	if !isempty(positional_names) && endswith(positional_names[end], "...")
		varargs_name = pop!(positional_names)[1:end-3]
		varargs_type = pop!(positional_types)
	end

	# A list of keyword argument names can be retrieved this way. Haven't
	# found a way to access keyword argument types yet.
	#kwarg_names = Base.kwarg_decl(method, typeof(methods.mt.kwsorter));

	return Subcommand(string(method.name), positional_names, positional_types,
		varargs_name, varargs_type)
end

function print_usage(cmd::Subcommand)
	print(stderr, "$(replace(cmd.name, '_' => ' '))")
	for arg in cmd.positional_names; print(stderr, " <$(arg)>"); end
	if cmd.varargs_name != ""; print(stderr, " <$(cmd.varargs_name)>..."); end
	#println(stderr, @doc(method.func).content[1].content[1].content[1])
	#if !isempty(kw_args)
	#	println(stderr, "\nOptions:")
	#end
end

function parse_option(value::AbstractString, typ::AbstractString)
	if startswith(typ, "Int") || startswith(typ, "Float")
		return parse(eval(Symbol(typ)), value)
	elseif typ == "Bool"
		lvalue = lowercase(value)
		if lvalue == "true" || lvalue == "yes"; return true; end
		if lvalue == "false" || lvalue == "no"; return false; end
	elseif typ == "IO"
		path = expanduser(value)
		if path == "-"; return stdin; end
		# TODO: Use GzipDecompressionStream instead?
		return endswith(path, ".gz") ? open(`gunzip -c $path`) : open(path)
	else
		return value
	end
end

# FIXME: This is a horrible hack, but necessary because Julia does
# not currently support introspection of keyword argument types.
function parse_speculatively(value::AbstractString)
	if lowercase(value) == "true"; return true; end
	if lowercase(value) == "false"; return false; end
	if all(isdigit, value); return parse(Int, value); end
	return try parse(Float64, value) catch; value; end
end

function execute_subcommand(cmd::Subcommand)
	name_parts = split(cmd.name, '_')
	if ARGS[1:min(length(name_parts), length(ARGS))] != name_parts
		return false
	end

	positional = []; varargs = []
	options = Dict{Symbol, Any}()

	for k in length(name_parts)+1:length(ARGS)
		if startswith(ARGS[k], "--")
			eq = findfirst(ARGS[k], '=')
			option = eq == nothing ? ARGS[k][3:end] : ARGS[k][3:eq-1]
			option = Symbol(replace(option, '-' => '_'))
			options[option] = eq == nothing ? true :
				parse_speculatively(ARGS[k][eq+1:end])
		else
			push!(positional, ARGS[k])
		end
	end

	if length(positional) < length(cmd.positional_names)
		print(stderr, "Usage: "); print_usage(cmd); print(stderr, "\n")
		return true
	end

	varargs = positional[length(cmd.positional_names)+1:end]
	positional = positional[1:length(cmd.positional_names)]

	if cmd.varargs_name != "" && length(varargs) == 0
		print(stderr, "Usage: "); print_usage(cmd); print(stderr, "\n")
		return true
	end

	positional = map(p -> parse_option(positional[p], cmd.positional_types[p]),
		1:length(positional))
	varargs = map(v -> parse_option(v, cmd.varargs_type), varargs)
	
	# FIXME: This assumes that subcommand handler is defined in Main
	eval(Meta.parse("Main.$(cmd.name)"))(positional..., varargs...; options...)

	# This does not work, Julia runtime claims that function does not support
	# keyword arguments (even if it does).
	#cmd.func(positional...; options...)
	return true
end

function subcommands(functions...)
	cmds = map(f -> Subcommand(f), functions)
	for cmd in cmds
		if execute_subcommand(cmd); return; end
	end

	println(stderr, "Available subcommands:")
	for cmd in cmds
		print(stderr, "  "); print_usage(cmd); print(stderr, "\n")
	end
end

end
