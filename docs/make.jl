using Pkg
using Documenter, DocumenterCitations, DemoCards
using Literate
using ParticleInCell2

@info "Generating examples using DemoCards"
examples_page, postprocess_democard_cb, demo_assets = makedemos("../examples")

@info "Gathering information from Project.toml"
PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ")
GITHUB = "https://github.com/adamslc/ParticleInCell2.jl"

@info "Making docs"
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"), style = :authoryear)

assets = ["assets/citations.css"]
isnothing(demo_assets) || (push!(assets, demo_assets))

makedocs(
    bib,
    authors = AUTHORS,
    sitename = "ParticleInCell2.jl Documentation",
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
    strict = get(ENV, "CI", nothing) == "true",
)

@info "Postprocessing DemoCards"
postprocess_democard_cb()

@info "Deploying docs"
deploydocs(repo = "github.com/adamslc/ParticleInCell2.jl.git")
