# transcription-model
Julia code used in J Rodriguez, G Ren, CR Day, K Zhao, CC Chow, and DR Larson. Intrinsic dynamics of an endogenous human gene reveal the basis of expression heterogeneity. Cell (in press)

fitscript.jl is a script that reads in data and fits model using a MCMC algorithm. Model data in the form of ON and OFF time distributions and steady state mRNA smFISH histograms is simulated with Gillespie algorithm or computed from direct solutions of the system Master equation
