# Stochastic "Wow"

Code and output files from 'Could the “Wow” signal have originated from a stochastic repeating beacon?' by Kipping & Gray (2022).

wow.f90 is a Fortran 90 program that changes the signal duration, T, and mean rate of signal recurrence, lambda. It slides through a grid of possible values, with indices k1 and k2, and outputs the likelihood of such a scenario. These results are saved and compiled in output_grid.dat. That file has some extra columns, and the columns are explained in the header.

mcmc.nb is a Mathematica notebook that runs a simple MCMC fit by generating a likelihood emulator from the likelihood grid and then walking in T and lambda space. You cna find the figures in the paper being generated inside that notebook. The output from three different MCMC runs are also provided in the repo: 1) mcmc_chain_thinner_OSU.dat 2) mcmc_chain_thinner_OSUHobart.dat 3) mcmc_chain_thinner_ALL.dat. These are essentially the same but we add additional constraints to the likelihood (as described in the paper) to include Hobart data in 2) and all available constraints in 3).
