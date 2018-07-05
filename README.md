# Pharmacokinetic-Pharmacodynamic Neuromuscular Junction Model Code

This repository is code for part of my dissertation project.  Electric field stimulation (EFS) is applied to tracheal ring segments. Force of contraction is measured experimentally and then computationally modeled as a series of pharmacokinetic compartments.

The raw data are contained in two folders corresponding to the upper trachea and the lower trachea. Both folders are subdivided into control (con) and capsaicin (cap) folders.  Data are CSV files formatted with "Time" in the far left hand column followed by increasing frequency response data from "0.1Hz" to "30Hz". The data are trachea responses measured in mg of tension. Time is in seconds with measurements every 0.02 seconds.  Total record time is 60 seconds.  All data have been formatted so that EFS begins at Time = 0. Data in the capsaicin folder has been pre-treated for 2 hours with 100 μM capsaicin and 10 μM indomethacin followed by a 30-minute washout period.  Not all data are created equal: in the inferior trachea only files 1, 2, 5, and 7 had characteristic responses and were used in analysis. Data files 1, 2, 3, 4, 5 were all used in the superior trachea analysis.

#### do.R and do_sup.R

The file do.R imports the necessary functions and the appropriate data and runs the model. Running this file will duplicate the approach I used to find the best fit parameters for the control and capsaicin inferior trachea data.  It will produce the output files that can be used in conjunction with the consensus_plots.R file to make plots and to find the aggregate parameters accross frequencies and tissues.

First it establishes a naive set of initial parameters and defines the First-Order model. Next it solves the model 100 times for the 0.1 Hz frequency and then for the 0.3, 1, 3, 10 Hz frequencies for both the capsaicin and control data. The 0.1 Hz frequencies are done separately because the experimental stimulation frequency is slighly different from 0.1 Hz and this difference can dramtically affect the fit.  As such the model fits the Frequency at 0.1 Hz.

The 0.1 Hz best fit is also used for the k~a~, k~e~, and pDose values for the complex models.  This model is solved next for both capsaicin and control. The final parameters for each tissue are averaged to determine a single set of best fit parameters for each tissue and treatment. Finally, it creates a second file with the results broken down by frequency rather than by tissue.

The file do_sup.R follows the same code, but performs it for the superior rather than the inferior trachea.

#### model_2drugs.R

This is a working file that contains instances of running the model functions described below. The model workflow begins with loading required packages and files, then the initial model parameters are defined, the selected pharmacokinetic model is compiled, and the data file to be tested is imported.  The model is run to find best-fit parameters and the results can be exported or graphed.

#### normalizedDF.R

The function `loadNormalizedDF` imports a selected data frame and normalizes the mg of tension response to a maximum. Parameters passed to the function include the section of trachea (lower = TRUE or FALSE) where the data was collected, data to be normalized ("con" or "cap" data), data to be used in the denominator to normalize (again, "con" or "cap" data), and the index for the file in its respective folder (integer that indicates the position of the file in the folder). In the lower trachea, capsaicin data is used to normalize both capsaicin and control data.  In the upper trachea however, control responses may have greater maxima than capsaicin responses and it may be necessary to normalize control responses to themselves.

#### defineModel.R

The `defineModel` function from defineModel.R creates and compiles the model to be tested using the `RxODE` package.  The constrictive ACh component of the model can be the First Order model ("simple") with first-order one-compartment  kinetics or our Inhibition model ("complex") with saturable kinetics for absorption and elimination. It is possible to create multiple neurotransmitter models with this function. The unknown relaxant component can be non-existent ("none") or present as a one-compartment first-order neurotransmitter ("simple"). If there is no unknown relaxant the effect model will default to a single neurotransmitter effect (pharmacodynamic) model.  If there is a defined unknown relaxant model there are two options for the way these two neurotransmitter effects are summed together. In the case of the "twoNT_1" option the relaxant effect is subtracted from the ACh effect and the product of the two effects is added back, analogous to the Bliss Independence model. In the "twoNT_2" option the unknown effect is instead subtracted from the ACh maximum effect. The function returns an object that contains the model, the models initial conditions, and the number of neurotransmitters in the model.  The information from this object is read by subsequent functions.

This function allows for the rapid selection from 6 relevent model choices using one line of code. Of course other models can be used, but these will have to be developed individually and have not been validated to work with the rest of the code presented here.

#### runModelFunctions.R

This file contains the major workhorses that bring the magic of the model together, the `run_mod1` and `final_drug_params` functions, as well as several other useful functions. After importing the data, defining a model, and choosing a set of initial parameters, `run_mod1` runs the model using the `RxODE` package.

`run_mod1` first maps the "init_params" to the parameters defined in the model. Next it builds the event table.  The event table is a list of doses of drug to be administered based on the EFS frequency. Finally it solves the model and returns a data frame of the compartment concentrations and predicted effects at each time point (every 0.02 s).

One run of `final_drug_params` finds a set of best-fit parameters for the chosen model.  The function runs `run_mod1` repeatedly, each time choosing new parameters and trying to minimize the sum of squares.  The parameters to be optimized are defined in the "bestfit" vector. Any combination of parameters from the "init_params" list can be fit by the function as long as they are in the defined model.  The results of `final_drug_params` can return either a data frame of either all the tested parameters (chosen = FALSE) or a data frame of only the tested parameters that were accepted as improved fits for the model (chosen = TRUE). In either case, the final row of the returned data frame are the best-fit parameters for that run.

The default number of iterations to find the best fit ("m") is 50 for rapid testing purposes, but 500 is more appropriate to ensure the model converges.  Two hyper-parameters ("lambda" and "testexp") can also adjust the rate of convergence.  The width (standard deviation) of the Gaussian curves used to sample new best-fit parameters is proportional to lambda.  A larger "lambda" will reduce the influence of the initial parameters by searching a larger parameter space, but may also waste many iterations by attempting numbers that have little chance of improving the fit.  The "testexp" hyper-parameter changes the penalty for having large variations from the experimental data.  The larger "testexp" is the more stringent the best-fit criteria are, but should always be an even integer (2 or 4 are recommended).  This function also imposes some hard boundaries on a few of the parameters.  AChE~max~ and M2~max~ and not allowed to go above 1.  The parameters pIC50~AChE~ and pIC50~M2~ are constrained to less than 10.

After the `final_drug_params` has been run `run_mod1` is run again on the final set of parameters to generate concentration vs. time data to be plotted.  An example of plotting the resulting data exists in model_2drugs.R using the  `ggplot2` R package. The example plotting function given plots the experimental data in black, initial results of `run_mod1` in violet, and the final best-fit results in blue.

`Iteration` is a wrapper function that takes a list of file indices (con_list) and some normalization inputs and generates n sets of final best-fit parameters for each file at each frequency (freq_list) using the defined model (ITmodel).   The default value of n used in most runs was 100. After each frequency it appends its results to a single CSV file per index. These files can be reimported later to calculate their statistics.

`facetgraph` is another wrapper function that produces a single multi-panel graph using `ggplot2` where each panel is a different frequency for the same model parameters. If consensus = FALSE, the init_params are taken as a starting point and `final_drug_params` is run a single time at each frequency . The resulting plot displays the experimental data in black, initial parameter results in blue, and the final best-fit results in red. Alternatively if consensus = TRUE, it will instead assume the "init_params" are the final consensus values and not run `final_drug_params`. Users can also supply consensus_params to plot a second set of values in the facet graph.  This is good for producing result graphs or for visual comparison of a model's fit at all frequencies.

The `aggregate_stats` function imports the results from `Iteration` and finds the mean, median and sd for each file given in the index by frequency. It exports them to a new file given by the output_filename.

The final function in runModelFunctions.R is `find_allfreq_params`.  This was developed after most of the data for the dissertation was generated.  Instead of finding best fit values for a single frequency like `final_drug_params`, it finds best fit values for all frequencies at once. The objective function it tries to minimize is the sum of the sum of squares/AUC. This still needs to be tested more rigorously, but could be a way to improve the fits of the model.

#### consesusPlot.R

consensusPlot.R is similar to model_drugs.R because it is a rough staging ground for arunning a collection of the functions discusses above. It uses the  `facetgraph` function, as well as finds the AUC, sums of squares (SS) and max results. It takes finalparameters.csv as input and can create plots usi. Each AUC, SS, and max data frame needs to be found manually and appended to the appropriate data frame. After they have all been run, then the results can be exported. This code is not completely organized and should be used with caution.
