# transcription-model
Julia code used in Rodriguez, Joseph, Gang Ren, Christopher R. Day, Keji Zhao, Carson C. Chow, and Daniel R. Larson. “Intrinsic Dynamics of a Human Gene Reveal the Basis of Expression Heterogeneity.” Cell 176, no. 1–2 (2019): 213-226.e18. https://doi.org/10.1016/j.cell.2018.11.026.  The code was originally written in Julia 0.5 but it has been updated to Julia 1.4.


fitscript.jl is a script that reads in data and fits model using a MCMC algorithm. Model data in the form of ON and OFF time distributions and steady state mRNA smFISH histograms is simulated with Gillespie algorithm or computed from direct solutions of the system Master equation.
