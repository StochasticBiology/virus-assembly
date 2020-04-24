# virus-assembly

Simulation code for virus capsid assembly

From Johnston et al. "Modelling the self-assembly of virus capsids", J Phys Condens Matt 22 104101 (2010) (free preprint here https://arxiv.org/pdf/0910.1916.pdf )

A Metropolis Monte Carlo scheme is used to simulate the thermodynamic interactions of different agents in a box with periodic boundary conditions. Agents can be pentamers, hexamers, or crowding agents, which interact through different potentials.

Command-line parameters can change the default parameterisation. By default, the size of bound clusters in the simulation and a current snapshot is output over time; time-labelled snapshots an also be output to allow animations to be produced.

This code is an odd C/C++ hybrid. Agents are objects with member functions, but no other object-orientated approaches are used. Compile with g++ (with maths library). This is a refactoring and amalgamation of several bits of old code and as such is a work in progress!

virus-sim.c -- wrapper code to be run, deals with I/O (command line arguments and file output) and implements the Monte Carlo simulation
virus-all.c -- describes agents as objects, providing member functions and other physical functions for the simulation
vectors-all.c -- provides overloaded operator and some functions for 3D vectors used in the simulation

run.sh -- bash script compiling the code and running a set of example physical experiments
get-stats.sh -- bash script summarising simulation outputs

virus-visualise.nb -- Mathematica notebook allowing visualisation of snapshots of simulation output
plot-expts.sh -- Gnuplot script plotting assembly yields with different parameters for the experiments in run.sh
