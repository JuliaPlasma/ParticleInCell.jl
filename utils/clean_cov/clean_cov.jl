directory = basename(pwd())
if !(directory in ["ParticleInCell2", "ParticleInCell2.jl"])
	@error "Script should be run the project root directory"
	exit(1)
end

using Pkg
Pkg.instantiate()

using Coverage
clean_folder("src")
clean_folder("test")
# clean_folder("ext")
