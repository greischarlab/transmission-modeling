Before running the code in any RMD files, open the R Project file ('Using models to identify the causes of 
pre-symptomatic transmission from human infection data.Rproj') at the root of the directory to correctly set the 
working directory.

# Literature review of CHI trials

'CHI trial lit review.RMD' contains the plotting code for Figure 1, which illustrates the data collected for 
our literature review of controlled human infection trials.

# Identifying correlates of pre-symptomatic transmission

'Fitting statistical model to human infection data.RMD' contains the code for all the cleaning and analysis of 
the data from Ge *et al*., 2023 and Zhou *et al*., 2023. This includes fitting the statistical model we use
to the individual shedding data, running linear regressions assessing the correlations between the statistical 
model parameters and the duration of pre-symptomatic transmission across all the participants, and bootstrapping 
the predicted values of the linear regressions to assess statistical uncertainty.

# Mechanistic within-host model

'Mechanistic model code.RMD' contains the code for the simulations and analyses of the within-host model we 
employ to produce a plausible scenario of pre-symptomatic transmission informed by earlier findings vis-a-vis 
data from the Ge *et al*., 2023 and Zhou *et al*., 2023 studies. 
