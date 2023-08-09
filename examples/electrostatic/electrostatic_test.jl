@testitem "Two-stream example" tags = [:integration, :example] begin
    include("two_stream.jl")
end

@testitem "beam-plasma example" tags = [:integration, :example] begin
    include("beam_plasma.jl")
end

@testitem "beam-cyclotron example" tags = [:integration, :example] begin
    include("beam_cyclotron.jl")
end

@testitem "Landau damping example" tags = [:integration, :example] begin
    include("landau_damping.jl")
end
