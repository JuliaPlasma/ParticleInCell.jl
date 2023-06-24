# Script modified from
# https://gist.github.com/pfitzseb/22493b0214276d3b65833232aa94bf11

using LanguageServer, StaticLint, SymbolServer

path = pwd()
root_file = path * "/" * "src/ParticleInCell2.jl"

s = LanguageServerInstance(Pipe(), stdout, path)
_, symbols = SymbolServer.getstore(s.symbol_server, path)
s.global_env.symbols = symbols
s.global_env.extended_methods = SymbolServer.collect_extended_methods(s.global_env.symbols)
s.global_env.project_deps = collect(keys(s.global_env.symbols))

f = StaticLint.loadfile(s, root_file)
StaticLint.semantic_pass(LanguageServer.getroot(f))

errors = Tuple{String,LanguageServer.Diagnostic}[]
for doc in LanguageServer.getdocuments_value(s)
    StaticLint.check_all(
        LanguageServer.getcst(doc),
        s.lint_options,
        LanguageServer.getenv(doc, s),
    )
    LanguageServer.mark_errors(doc, doc.diagnostics)

    for d in doc.diagnostics
        push!(errors, (doc._path, d))
    end
end

for (file, error) in errors
    @warn error.code msg = (error.message) _file = file _line = (error.range.start.line)
end
