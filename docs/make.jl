using Pkg, Documenter, DocumenterCitations, Literate, ParticleInCell2

@info "Generating tutorial materials from Literate script"
tutorial_path = joinpath(@__DIR__, "literate_src", "tutorial.jl")
output_path = joinpath(@__DIR__, "src")
Literate.markdown(tutorial_path, output_path, documenter = true)
Literate.notebook(tutorial_path, output_path, documenter = true)
Literate.script(tutorial_path, output_path, documenter = true)

@info "Gathering information from Project.toml"
PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ")
GITHUB = "https://github.com/adamslc/ParticleInCell2.jl"

@info "Making docs"
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"), style = :authoryear)

makedocs(
    bib,
    authors = AUTHORS,
    sitename = "ParticleInCell2.jl Documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = String["assets/citations.css"],
        footer = "[$NAME.jl]($GITHUB) v$VERSION",
    ),
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "examples/index.md",
        "Plasma simulation theory" => [
            "Introduction" => "theory/index.md",
            "PIC simulation" => "theory/intro_to_pic.md",
        ],
        "reference/index.md",
        "references.md",
    ],
    strict = get(ENV, "CI", nothing) == "true",
)

@info "Deploying docs"
deploydocs(repo = "github.com/adamslc/ParticleInCell2.jl.git")
