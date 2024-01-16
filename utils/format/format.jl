directory = basename(pwd())
if !(directory in ["ParticleInCell", "ParticleInCell.jl"])
    @error "Format script should be run the project root directory"
    exit(1)
end

using Pkg
Pkg.instantiate()

using JuliaFormatter
format(".", verbose = true)
