# Monte Carlo Coalescence model

## General structure of the project

The coalescence model is basically an order of 10 numerical integral, depending on whether it is for mesons or barions. These integrals are performed via Monte Carlo integration. This integration is called through the function called *rapidity_momentum_distribution( func_vars_array , par_array , integ_lims_array , N_iter )* which has 3 argument at the moment. The 2 vector arguments are the functions variables, i.e., pT and y in our case, the parameters of the function, like the mass or the temperature. The third argument is a **long int** that pass the number of iteration to the Monte Carlo algorithm. It should be given as a command line argument.

The detailed structure of the model can be studied in the [overleaf document](https://www.overleaf.com/project/622f15b1b8ea1e1d0f2566c8) that I maintain. In general, the function that is to be integrated consists two main parts: the coalescence distribution and the momentum distributions of the partons. The latter could be two or three separate distribution with a similar structure depending on the studied case, i.e., mesons or barions are modelled. The coalescence distribution is a probability-like distribution in the phase space. It describes how probable is that the two/three considered parton coalesce into a meson/barion. The momentum dsitributions are taken as thermal distributions below p_T < 2 GeV/c and power-law above that.

The data folder contains experimental data from various sources that are listed in the subsequent [README file]().