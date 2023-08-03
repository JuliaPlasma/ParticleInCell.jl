using Pkg
using Documenter, DocumenterCitations, DemoCards
using Literate
using ParticleInCell2

@info "Generating tutorial materials from Literate script"
tutorial_path = joinpath(@__DIR__, "literate_src", "tutorial.jl")
tutorial_output_path = joinpath(@__DIR__, "src")
Literate.markdown(tutorial_path, tutorial_output_path, documenter = true)
Literate.notebook(tutorial_path, tutorial_output_path, documenter = true)
Literate.script(tutorial_path, tutorial_output_path, documenter = true)

@info "Generating examples using DemoCards"
examples_page, postprocess_democard_cb, demo_assets = makedemos("literate_src/examples")

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
        "Tutorial" => "tutorial.md",
        examples_page,
        "Plasma simulation theory" => [
            "Introduction" => "theory/index.md",
            "PIC simulation" => "theory/intro_to_pic.md",
        ],
        "reference/index.md",
        "references.md",
    ],
    strict = get(ENV, "CI", nothing) == "true",
)

@info "Postprocessing DemoCards"
postprocess_democard_cb()

@info "Clean up generated tutorial materials"
rm(jointpath(tutorial_output_path, "tutorial.md"))
rm(jointpath(tutorial_output_path, "tutorial.jl"))
rm(jointpath(tutorial_output_path, "tutorial.ipynb"))

@info "Deploying docs"
deploydocs(repo = "github.com/adamslc/ParticleInCell2.jl.git")
