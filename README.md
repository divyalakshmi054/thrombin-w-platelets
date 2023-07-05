# BST Model for Thrombin Generation in Whole Blood

To run simulations: Assign visit number in `sample_ensemble.jl`, then `include("sample_ensemble.jl")`

To estimate parameters by reducing differences: Assign parameter guesses, then `include("estimate_model_parameters.jl")`

To use estimated parameters to generate other cases: Pick a PSET in `construction.jl`. Then,  `include("construction.jl")`. This file has a plotting routine to generate FIIa v. Time *for each patient* and store to `figs` in pwd.

To generate average plots for N patients: Open the appropriate experimental data file and in a loop running from 1 to N, open simulation file(s) in `plot_thrombin.jl` to compute mean and standard error. Then, `include("plot_thrombin.jl")`

To run Sobol sensitivity analysis: Edit number of samples and bootreps in gsa method in `sensitivity_analysis.jl`, then `include("sensitivity.jl")`. This routine also dumps Sobol results to disk to a `sobol` directory within pwd

To plot Sobol indices: `include("plot_sensitivity.jl")`; you may need to update the CSV filename or the path to file containing the sobol results
