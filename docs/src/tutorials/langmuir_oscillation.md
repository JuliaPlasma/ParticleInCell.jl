# Introduction to PIC methods: Langmuir oscillations
In this first tutorial, we introduce the basic concepts of PIC simulation, and
then demonstrate how PIC simulation can be used to model an electrostatic
Langmuir oscillation in one dimension. If you are already familar with PIC
methods, you can skip to [the code](TODO), to see how the simulation is
implemented in `ParticleInCell2.jl`.

## Particle-in-Cell basics
Classical plasma physics considers the motion of charged particles. The
dynamics of these particles will be effected by the presence of externally
imposed electric and magnetic fields, which is relatively easy to model since
the motion of each particle is independent. However, the particles will
themselves source an electric field---and, if they are moving fast enough, a
magnetic field---due to their charge. The dynamics of other particles will be
influenced by these 'self-consistent' fields, which corresponding source fields
of their own. Thus, the dynamics of all of the particles is coupled, and so
their equations of motion must be solve together. For a typical plasma, the
resulting differential equations cannot be solved analytically, and the
number of degrees of freedom means that direct computational integration of
the equations is also impossible. Thus, plasma physics relies on various
simplifications and assumptions to create reduced models that can be solved---
either exactly, or approximately.

### Approximations to the true particle distribution function
It is convenient to represent the locations and momenta of the particles using
the distribution function
```math
f(x,p,t) = \sum_{i=1}^N \delta(x - x_i(t)) \delta(p - p_i(t)),
```
where $N$ is the total number of particles. The evolution of this distribution
function can be written using the Maxwell-Boltzmann system
```math
TODO.
TODO-Maxwell
```
Unfortunately, for almost all plasmas of interest, $N$ is enormous; thus direct
simulation of the equations of motion is computationally intractable. The
solution is to recognize that if $N$ is large, then any physically small
phase-space region will contain many particles, and so we c

#### Todo
- course-graining
- collisionless plasmas
- briefly mention thermalization, two-fluid, and MHD

## Discretized solutions of the Boltzmann-Poisson equations
We have seen that a kinetic plasma can be approximated using the
Vlasov-Boltzmann equation, along with an appropriate equation of motion for the
electromagnetic fields. However, the resulting equation is still not easy to
analyse for an arbitrary $f$. We therefore turn to computational simulation of
the equations of motion to understand the plasma dynamics.

In order to simulate the plasma, we must first discretize the equations so they
can be represented on a computer. One option, called a Vlasov code, represents
the distribution function as

## The standard PIC cycle

## Analytical derivation of Langmuir oscillations

## Simulation of Langmuir oscillations

