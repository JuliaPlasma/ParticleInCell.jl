# Plasma Simulation Theory

The fundamental physics governing the dynamics of a plasma have been well
understood for over a century, and yet plasma physics remains an active area of
research. This is because the dynamics of a plasma are highly nonlinear, and
it is therefore difficult to make analytic statements about how a given plasma
will behave over long time-spans. Instead, theoretical plasma physicists often
rely on simulation to understand plasma dynamics, and to make general
statements about the behaviors of particular classes of plasmas.

Unfortunately, the simulation of plasmas is itself a quite challenging task
because plasmas are composed of *many*, *many*, charged particles. As a result,
there are far too many degrees of freedom to exactly solve the full equations
of motion for a given plasma. Instead, physicists rely on approximations to
derive physically relevant models for the systems in question.

The most "realistic" class of models are called *kinetic* models. In these
models, the individual charged particles of each species are averaged to create
distribution functions that depends on both position and velocity. The
distribution functions give the likelihood of finding a charged particle of a
particular species at any point in phase space. These distribution functions,
along with a method for calculating the self-consistent interaction between the
particles, yields a set of approximate equations of motion describing the
plasma dynamics.

If the particle species are nonrelativistic, then the particles are not
influenced by self-consistent magnetic fields, and so the interactions can be
modeled using electrostatics. However, once the particles (typically the
lightweight electrons) become relativistic, the sourced magnetic field must be
taken into account, and so full electromagnetic interactions are required.

Over long periods of time, the plasma will begin to thermalize---the
distribution of particle velocities will become closer and closer to
Maxwellian. This fact can be used to drastically improve the size and speed of
simulations by using the two-fluid and magnetohydrodynamic (MHD) approximations
of plasma dynamics. As this package does not implement algorithms for these
simulation methods, we will not discuss the details of these methods further.

In the following sections, we describe the details of particle-in-cell
algorithms, a specific type of kinetic simulation algorithm. For a more in
depth introduction of PIC simulation theory, see the books by
[birdsall2004](@citet) and [hockney1989](@citet).
