using Pkg
using Documenter, DocumenterCitations, DemoCards
using Literate
using ParticleInCell

@info "Generating examples using DemoCards"
examples_page, postprocess_democard_cb, demo_assets =
    makedemos("../examples", filter_function = x -> !endswith(x.path, "_test.jl"))

@info "Gathering information from Project.toml"
PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ")
GITHUB = "https://github.com/JuliaPlasma/ParticleInCell.jl"

@info "Making docs"
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"), style = :authoryear)

assets = ["assets/citations.css"]
isnothing(demo_assets) || (push!(assets, demo_assets))

makedocs(
    plugins = [bib],
    authors = AUTHORS,
    sitename = "ParticleInCell.jl Documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = assets,
        footer = "[$NAME.jl]($GITHUB) v$VERSION",
    ),
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => "examples/tutorial/langmuir_oscillation.md",
        examples_page,
        "Plasma simulation theory" => [
            "Introduction" => "theory/index.md",
            "PIC simulation" => "theory/intro_to_pic.md",
        ],
        "manual/index.md",
        "references.md",
    ],
    modules = [ParticleInCell],
)

@info "Postprocessing DemoCards"
postprocess_democard_cb()

@info "Deploying docs"
deploydocs(repo = "github.com/JuliaPlasma/ParticleInCell.jl.git")
