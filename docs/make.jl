using Documenter, ParticleInCell2

makedocs(
    sitename = "ParticleInCell2.jl Documentation",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "examples/index.md",
        "theory/index.md",
        "reference/index.md",
    ],
)

deploydocs(repo = "github.com/adamslc/ParticleInCell2.jl.git")
