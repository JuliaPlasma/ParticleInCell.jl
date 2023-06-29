@testitem "Aqua" tags = [:basics] begin
    using Aqua
    # The unbound_args test currently has a false positive, because it thinks the grid
    # might be zero dimensional, in wich case T and U would be unbound. This could be
    # fixed by using a construction like Tuple{T, Vararg{T, D1}} where D1 = D-1, but
    # this is annoying becasue I want the value of D as a type parameter in several
    # places. For now, I have just disabled the test.
    Aqua.test_all(ParticleInCell2, unbound_args = false)
end
