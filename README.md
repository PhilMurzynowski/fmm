Files for the Fast Multipole Method.

Simulate and visualize by running `Simulate.jl`. Number of bodies, precision, size of time steps, number of time steps, gravitational constant, etc. can all be changed within `Simulate.jl`. Once all necessary packages are installed should only take at most a few seconds to boot the first run, and then the simulation is set to run for approximately 10 seconds.

Quadtree.jl stores the the quadtree struct and associated functions.
Multipole.jl contains all the functions necessary to manipulate expansions.
Simulate.jl can be run for a real time n-body simulation using a configurable about of particles, expansion terms, gravity constant, etc.
Accuracy.jl checks the accuracy of the fast multipole method for a given set of particles across a range of possible expansion sizes.
FlopCount.jl counts the flops in the FMM across a range of number of expansion terms and number of bodies.
Naive.jl is the direct O(N^2) implementation.
Timing.jl times execution for different numbers of expansion coefficients and potentially number of bodies.
TimingExperiments.jl holds some small scale timing experiments which larger optimizations were later based on (power operation, etc.).
BinomialTable.jl creates binomial lookup tables with different  access patterns.
FourColorSort.jl implements a four color version of the Dutch national flag algorithm along with variations for convenience.

The barnes hut directory contains a port of ParallelBarnesHut.jl for 2D written serially to avoid any potential issues with comparing flops,
and similar accuracy, simulation, and flop counting setups.
