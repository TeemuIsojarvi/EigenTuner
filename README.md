# EigenTuner
A C++ code for partial solution of the quantum inverse scattering problem for a finite number of eigenvalues.

The program begins with a 1D finite square well potential energy, with length L and depth V0 given as an input, and iteratively adds potential energy perturbation terms that make the lowest N energy eigenvalues converge towards user-inputted target values. At least up to 6 energies can be effectively tuned to wanted locations with the current version.

This program has been tested in Ubuntu Linux and compiled as "g++ -O3 -ffast-math eigentuner.cpp" or "g++ -O3 -ffast-math eigentuner_sqr_density.cpp".

The code in its current form contains a lot of global variables and is not very optimized for performance, but this will be corrected in later versions.

More information about the application can be found from the documentation file "Documentation_v0.95.pdf". Some questions about this topic are also considered in my Wordpress blog: https://compphysdiary.wordpress.com/ .
