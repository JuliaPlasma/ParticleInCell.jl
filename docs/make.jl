using Documenter, Literate, ParticleInCell2

@info "Generating tutorial materials from Literate script"
tutorial_path = joinpath(@__DIR__, "literate_src", "tutorial.jl")
output_path = joinpath(@__DIR__, "src")
Literate.markdown(tutorial_path, output_path, documenter = true)
Literate.notebook(tutorial_path, output_path, documenter = true)
Literate.script(tutorial_path, output_path, documenter = true)

@info "Making docs"
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

@info "Deploying docs"
deploydocs(repo = "github.com/adamslc/ParticleInCell2.jl.git")
