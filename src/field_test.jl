@testitem "Field" tags = [:unit] begin
    grid = UniformCartesianGrid(
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (10, 10, 10),
        (true, true, true),
    )
    scalar_node_field = Field(grid, ParticleInCell.node, 1, 3)
    vector_node_field = Field(grid, ParticleInCell.node, 3)
    vector_edge_field = Field(grid, ParticleInCell.edge, 3)
    vector_face_field = Field(grid, ParticleInCell.face, 3)

    @test num_elements(scalar_node_field) == 1
    @test num_elements(vector_node_field) == 3
    @test num_elements(vector_edge_field) == 3
    @test num_elements(vector_face_field) == 3

    @test eachindex(scalar_node_field) == CartesianIndices((4:13, 4:13, 4:13))
    @test eachindex(vector_node_field) == CartesianIndices((1:10, 1:10, 1:10))

    @test ParticleInCell.cell_index_to_cell_coords(scalar_node_field, CartesianIndex(4, 4, 4)) == (0, 0, 0)
    @test ParticleInCell.cell_index_to_cell_coords(scalar_node_field, CartesianIndex(14, 14, 14)) == (10, 10, 10)
    @test ParticleInCell.cell_index_to_cell_coords(scalar_node_field, CartesianIndex(4, 4, 4), 1) == (0, 0, 0)
    @test ParticleInCell.cell_index_to_cell_coords(scalar_node_field, CartesianIndex(14, 14, 14), 3) == (10, 10, 10)

    @test ParticleInCell.cell_index_to_cell_coords(vector_node_field, CartesianIndex(1, 1, 1)) == (0, 0, 0)
    @test ParticleInCell.cell_index_to_cell_coords(vector_node_field, CartesianIndex(11, 11, 11)) == (10, 10, 10)
    @test ParticleInCell.cell_index_to_cell_coords(vector_node_field, CartesianIndex(1, 1, 1), 1) == (0, 0, 0)
    @test ParticleInCell.cell_index_to_cell_coords(vector_node_field, CartesianIndex(11, 11, 11), 3) == (10, 10, 10)

    @test ParticleInCell.cell_index_to_cell_coords(vector_edge_field, CartesianIndex(1, 1, 1)) == (0, 0, 0)
    @test ParticleInCell.cell_index_to_cell_coords(vector_edge_field, CartesianIndex(10, 10, 10)) == (9, 9, 9)
    @test ParticleInCell.cell_index_to_cell_coords(vector_edge_field, CartesianIndex(1, 1, 1), 1) == (0.5, 0, 0)
    @test ParticleInCell.cell_index_to_cell_coords(vector_edge_field, CartesianIndex(10, 10, 10), 2) == (9, 9.5, 9)
    @test ParticleInCell.cell_index_to_cell_coords(vector_edge_field, CartesianIndex(11, 11, 11), 3) == (10, 10, 10.5)

    @test ParticleInCell.cell_index_to_cell_coords(vector_face_field, CartesianIndex(1, 1, 1)) == (0, 0, 0)
    @test ParticleInCell.cell_index_to_cell_coords(vector_face_field, CartesianIndex(10, 10, 10)) == (9, 9, 9)
    @test ParticleInCell.cell_index_to_cell_coords(vector_face_field, CartesianIndex(1, 1, 1), 1) == (0, 0.5, 0.5)
    @test ParticleInCell.cell_index_to_cell_coords(vector_face_field, CartesianIndex(10, 10, 10), 2) == (9.5, 9, 9.5)
    @test ParticleInCell.cell_index_to_cell_coords(vector_face_field, CartesianIndex(11, 11, 11), 3) == (10.5, 10.5, 10)

    @test ParticleInCell.cell_coords_to_cell_index(scalar_node_field, (0, 0, 0)) == CartesianIndex(4, 4, 4)
    @test ParticleInCell.cell_coords_to_cell_index(scalar_node_field, (10.0, 10.0, 10.0)) == CartesianIndex(14, 14, 14)

    @test ParticleInCell.cell_coords_to_cell_index(vector_edge_field, (0, 0, 0)) == CartesianIndex(1, 1, 1)
    @test ParticleInCell.cell_coords_to_cell_index(vector_edge_field, (10.0, 10.0, 10.0)) == CartesianIndex(11, 11, 11)
end
