
include("mcmc.jl")    # MCMC simulator fitting FISH and Dwell time distributions simultaneously
include("dwelltimeCoupled.jl")  # Closed form solution of ON and OFF time distributions for coupled allele system
include("GillespieSimulator.jl")  # Gillespie simulator for LC and smFISH distributions
include("chicorrburstcond.jl")  # Gillespie simulator for correlations

using Distributions
using StatsBase

samples = 10   # Number of mcmc steps
scheme = 1    # 0 uses closed form LC computation, 1 uses Gillespie simulator

nAlleles = 5  # Number for coupled alleles
nGstates = 3  # Number of G states
nRsteps = 2   # Number of R steps

# total number of params
nparams = 2*(nGstates -1) + nRsteps + 3

LCsets = [4;5;6;7]  # Indices of parameter columns corresponding Live Cell data sets
FISHsets = [1;2;3]  # Indices of parameter columns corresponding to smFISH data sets
doseset = ([1],[4],[5],[2;6],[3;7])  # Indices of shared doses, so parameters for those models are identical

fixedeffects = [1;2;4;5;6;7;9]  # Parameter indices are the same for all doses
randomeffects = [3] # Parameter indices that depend on dose

nbursts = 2000 # total number of bursts for LC simulations
total =  2000000  # maximum number of steps for LC Gillespie simulations
tmax = 1e8  #  total time duration for smFISH simulation
totalf = 2000000  # maximum number of steps for smFISH Gillespie simulation
nbchor = 400 # total number of runs for allele cross-correlation estimates

temp = [1.;1.] # Temperature for mcmc

r = readdlm("rateparamsG3R2.txt",',')  #Input initial condition for parameters

# Four Live Cell data sets
dataLC = Array{Array{Float64,2},1}(4)
dataLC[1] = readdlm("dataLCMEM.txt",',')
dataLC[2] = readdlm("dataLCdeletion.txt",',')
dataLC[3] = readdlm("dataLCdose.05.txt",',')
dataLC[4] = readdlm("dataLCdose.5.txt",',')

# Three smFISH data sets
dataFISH = Array{Array{Int,2},1}(3)
dataFISH[1] = readdlm("dataFISHdose0.txt")
dataFISH[2] = readdlm("dataFISHdose.05.txt")
dataFISH[3] = readdlm("dataFISHdose.5.txt")


runoutput = mcmc(dataLC,dataFISH,r,LCsets,FISHsets,fixedeffects,randomeffects,doseset,nGstates,nRsteps,nAlleles,scheme,samples,total,totalf,tmax,nbursts,nbchor,temp);
