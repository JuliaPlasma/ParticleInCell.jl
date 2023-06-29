var documenterSearchIndex = {"docs":
[{"location":"tutorials/langmuir_oscillation/#Introduction-to-PIC-methods:-Langmuir-oscillations","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"","category":"section"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"In this first tutorial, we introduce the basic concepts of PIC simulation, and then demonstrate how PIC simulation can be used to model an electrostatic Langmuir oscillation in one dimension. If you are already familar with PIC methods, you can skip to the code, to see how the simulation is implemented in ParticleInCell2.jl.","category":"page"},{"location":"tutorials/langmuir_oscillation/#Particle-in-Cell-basics","page":"Introduction to PIC methods: Langmuir oscillations","title":"Particle-in-Cell basics","text":"","category":"section"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"Classical plasma physics considers the motion of charged particles. The dynamics of these particles will be effected by the presence of externally imposed electric and magnetic fields, which is relatively easy to model since the motion of each particle is independent. However, the particles will themselves source an electric field–-and, if they are moving fast enough, a magnetic field–-due to their charge. The dynamics of other particles will be influenced by these 'self-consistent' fields, which corresponding source fields of their own. Thus, the dynamics of all of the particles is coupled, and so their equations of motion must be solve together. For a typical plasma, the resulting differential equations cannot be solved analytically, and the number of degrees of freedom means that direct computational integration of the equations is also impossible. Thus, plasma physics relies on various simplifications and assumptions to create reduced models that can be solved–- either exactly, or approximately.","category":"page"},{"location":"tutorials/langmuir_oscillation/#Approximations-to-the-true-particle-distribution-function","page":"Introduction to PIC methods: Langmuir oscillations","title":"Approximations to the true particle distribution function","text":"","category":"section"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"It is convenient to represent the locations and momenta of the particles using the distribution function","category":"page"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"f(xpt) = sum_i=1^N delta(x - x_i(t)) delta(p - p_i(t))","category":"page"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"where N is the total number of particles. The evolution of this distribution function can be written using the Maxwell-Boltzmann system","category":"page"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"TODO\nTODO-Maxwell","category":"page"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"Unfortunately, for almost all plasmas of interest, N is enormous; thus direct simulation of the equations of motion is computationally intractable. The solution is to recognize that if N is large, then any physically small phase-space region will contain many particles, and so we c","category":"page"},{"location":"tutorials/langmuir_oscillation/#Todo","page":"Introduction to PIC methods: Langmuir oscillations","title":"Todo","text":"","category":"section"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"course-graining\ncollisionless plasmas\nbriefly mention thermalization, two-fluid, and MHD","category":"page"},{"location":"tutorials/langmuir_oscillation/#Discretized-solutions-of-the-Boltzmann-Poisson-equations","page":"Introduction to PIC methods: Langmuir oscillations","title":"Discretized solutions of the Boltzmann-Poisson equations","text":"","category":"section"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"We have seen that a kinetic plasma can be approximated using the Vlasov-Boltzmann equation, along with an appropriate equation of motion for the electromagnetic fields. However, the resulting equation is still not easy to analyse for an arbitrary f. We therefore turn to computational simulation of the equations of motion to understand the plasma dynamics.","category":"page"},{"location":"tutorials/langmuir_oscillation/","page":"Introduction to PIC methods: Langmuir oscillations","title":"Introduction to PIC methods: Langmuir oscillations","text":"In order to simulate the plasma, we must first discretize the equations so they can be represented on a computer. One option, called a Vlasov code, represents the distribution function as","category":"page"},{"location":"tutorials/langmuir_oscillation/#The-standard-PIC-cycle","page":"Introduction to PIC methods: Langmuir oscillations","title":"The standard PIC cycle","text":"","category":"section"},{"location":"tutorials/langmuir_oscillation/#Analytical-derivation-of-Langmuir-oscillations","page":"Introduction to PIC methods: Langmuir oscillations","title":"Analytical derivation of Langmuir oscillations","text":"","category":"section"},{"location":"tutorials/langmuir_oscillation/#Simulation-of-Langmuir-oscillations","page":"Introduction to PIC methods: Langmuir oscillations","title":"Simulation of Langmuir oscillations","text":"","category":"section"},{"location":"reference/#Reference-Manual","page":"Reference Manual","title":"Reference Manual","text":"","category":"section"},{"location":"reference/","page":"Reference Manual","title":"Reference Manual","text":"Modules = [ParticleInCell2]","category":"page"},{"location":"reference/#ParticleInCell2.AbstractGrid","page":"Reference Manual","title":"ParticleInCell2.AbstractGrid","text":"Parent type for all grid objects, which are used to define the simulation domain, and to convert between coordinate systems. There are three different numbering systems that can refer to a location in the simulation domain:\n\nThe 'physical coordinates' of a point are the real (dimensionalful) coordinates associated with that point. This value can range from the lower bounds to the upper bounds of the simulation. This value will typically take the form Vector{T} or NTuple{N, T} where T <: Real.\nThe 'cell coordinates' of a point is the non-dimensional location of the point in units of cell lengths. This value can range from 0 to num_cells - 1, or outside this range if guard cells are included. The value will typically have the type NTuple{N, Int}.\nThe 'cell index' of a point is the CartesianIndex that can be used to index into field arrays at that point. This value must strictly be confined to axes(field.values), which, for any given dimension, will typically range from 1 to numcells + 2*numguard_cells + 1.\n\nThe first two types of indexing, physcoords and cellcoords, are independent of the number of guard cells in a given field, and depend only on grid quantities. Thus utilities for converting between these systems require only a reference to a grid object. On the other hand, the utilities for cell_index are specific the field being used, and thus those must be provided an AbstractField to do the coordinate conversion.\n\nIn general, physical coordinates are useful when considering the location of a particle, while the cell index is used to interpolate to and from the particle locations. The cell coordinates are useful for some interpolation algorithms, especially those that are defined for non-uniform grids.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ParticleInCell2.UniformCartesianGrid","page":"Reference Manual","title":"ParticleInCell2.UniformCartesianGrid","text":"UniformCartesianGrid(lower_bounds, upper_bounds, num_cells, periodic)\n\nThe simplest grid type, which represents a set of equally spaced rectangular cells. The grid can optionally be periodic in one or more dimensions.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ParticleInCell2.cell_coords_to_phys_coords","page":"Reference Manual","title":"ParticleInCell2.cell_coords_to_phys_coords","text":"cell_coords_to_phys_coords(grid, idxs, [offset, component])\n\nConverts the cell coordinates idxs to a physical coordinate using the geometry specified in grid. Optionally, an offset and component can be specified to get the physical coordinates of a specific edge or face of the cell. The component argument is one-indexed.\n\nFor more information on the different types of coordinate systems, see AbstractGrid.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ParticleInCell2.cell_lengths","page":"Reference Manual","title":"ParticleInCell2.cell_lengths","text":"cell_lengths(grid, [cell_coords])\n\nReturns the length of the cell located at cell_coords. For uniform grids, the cell_coords argument is optional.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ParticleInCell2.compute_bspline_coeffs","page":"Reference Manual","title":"ParticleInCell2.compute_bspline_coeffs","text":"compute_bspline_coeffs(degree, [T])\n\nReturns a vector of vectors of the coefficients for the b-spline with polynomial degree, and unit-spaced knots, centered at zero.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ParticleInCell2.compute_knots-Tuple{Any}","page":"Reference Manual","title":"ParticleInCell2.compute_knots","text":"compute_knots(degree)\n\nReturns a vector of the knot locations for a polynomial b-spline with degree.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ParticleInCell2.phys_coords_to_cell_coords","page":"Reference Manual","title":"ParticleInCell2.phys_coords_to_cell_coords","text":"phys_coords_to_cell_coords(grid, xs)\n\nConverts the physical coordinate xs to a grid coordinate using the geometry specified in grid.\n\nFor more information on the different types of coordinate systems, see AbstractGrid.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ParticleInCell2.sim_lengths","page":"Reference Manual","title":"ParticleInCell2.sim_lengths","text":"sim_lengths(grid)\n\nReturns a tuple of the length of the simulation domain in each dimension.\n\n\n\n\n\n","category":"function"},{"location":"tutorials/#Tutorials","page":"Listing of tutorials","title":"Tutorials","text":"","category":"section"},{"location":"tutorials/","page":"Listing of tutorials","title":"Listing of tutorials","text":"In this section, we provide a few tutorials to get you up and running with ParticleInCell2. Each tutorial walks step-by-step through creating and running a PIC simulation, and then using the results to measure some physical quantity.","category":"page"},{"location":"tutorials/#Listing-of-tutorials","page":"Listing of tutorials","title":"Listing of tutorials","text":"","category":"section"},{"location":"tutorials/","page":"Listing of tutorials","title":"Listing of tutorials","text":"Introduction to PIC methods: Langmuir oscillations","category":"page"},{"location":"#ParticleInCell2.jl:-Kinetic-Plasma-Physics-in-Julia","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","text":"","category":"section"},{"location":"","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","text":"ParticleInCell2 is a package for doing plasma physics in Julia. Specifically, it focuses on the simulation of kinetic (non-thermal) plasmas using particle-in-cell (PIC) algorithms.","category":"page"},{"location":"#Installation","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"Installation","text":"","category":"section"},{"location":"","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","text":"ParticleInCell2 is currently not registered in the Julia package registry. Thus, to install this package, you should dev it:","category":"page"},{"location":"","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","text":"using Pkg\nPkg.develop(url=\"https://github.com/adamslc/ParticleInCell2.jl\")","category":"page"},{"location":"","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","text":"Once installed, you can load the package using","category":"page"},{"location":"","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","text":"using ParticleInCell2","category":"page"},{"location":"#Documentation-sections","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"Documentation sections","text":"","category":"section"},{"location":"","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","text":"The ParticleInCell2 documentation is divided into the following sections:","category":"page"},{"location":"","page":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","title":"ParticleInCell2.jl: Kinetic Plasma Physics in Julia","text":"Tutorials: step-by-step guides for basic package functionality.\nHow-to Guides: more complex examples of package capabilities.\nDiscussions: overviews of plasma simulation in general, and the specific design of the ParticleInCell2 package.\nReference Manual: detailed description of each public function and type in the package.","category":"page"},{"location":"discussion/#Discussions","page":"Discussions","title":"Discussions","text":"","category":"section"},{"location":"how_to/#How-to-Guides","page":"How-to Guides","title":"How-to Guides","text":"","category":"section"},{"location":"how_to/","page":"How-to Guides","title":"How-to Guides","text":"This section contains instructions for running more complex simulations.","category":"page"}]
}