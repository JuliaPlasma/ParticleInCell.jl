directory = basename(pwd())
if !(directory in ["ParticleInCell2", "ParticleInCell2.jl"])
    @error "Format script should be run the project root directory"
    exit(1)
end

using Pkg
Pkg.instantiate()

using JuliaFormatter
format(".", verbose=true)
