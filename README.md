# ParticleInCell2

## Todo:
 * Add examples
   * Two-stream instability (check frequency)
 * Add LinearFieldSolve step
   * Figure out a flexible way to enable stencil definitions
 * Add dump/restore abilities. Possible interface
   * dump(obj, dump_prefix::String, dump_num::Int) must produce a file
         {dump_prefix}_{obj_name}_{dump_num}.h5
     If an object does not change state over the course of the simulation, it is
     ok to ignore further dumps.
   * restore(::Type{ObjType}, dump_prefix::String, dump_num::Int)
     restore(::Type{ObjType}, file::H5File)
     Objects should be required to write there type to an attribute in the H5
     file, so that a restore program can read it
   * How much should I require the simulation to be restored only from on disk
     data vs require the original simulation file to be avaliable too?
 * AbstractUpdateStep interface:
   * step!(obj, dt) -- should I include dt? Should definitely get rid of sim
   * dump
   * restore
   * Need to come up with some nice way to systematically test these...
 * Add the ability to have subdomains
   * Locally, subdomains should be processed in parallel using Julia's coroutine,
   * Should also have 'ghost subdomains' that refer to a subdomain that is not on
     the current node. Then the communicate steps can use MPI to communicate ghost
     cells off node.
