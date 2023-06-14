#!/bin/sh

julia --project=utils/format -e '
directory = basename(pwd())
if !(directory in ["ParticleInCell2", "ParticleInCell2.jl"])
	@error "Format script should be run the project root directory"
	exit(1)
end

using Pkg
Pkg.instantiate()
Pkg.activate(".")
include("utils/lint/lint.")'
