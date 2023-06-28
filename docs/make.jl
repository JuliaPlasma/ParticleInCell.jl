using Documenter, ParticleInCell2

makedocs(
    sitename = "ParticleInCell2.jl Documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "index.md",
        "Tutorials" => [
            "Listing of tutorials" => "tutorials/index.md",
            "tutorials/langmuir_oscillation.md"
        ],
        "how_to/index.md",
        "discussion/index.md",
        "reference/index.md",
    ],
)

deploydocs(
    repo = "github.com/adamslc/ParticleInCell2.jl.git",
)
