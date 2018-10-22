
function mcmc(dataLC::Array{Array{Float64,2},1},dataFISH::Array{Array{Int,2},1},rin::Array{Float64,2},LCsets::Array{Int,1},FISHsets::Array{Int,1},fixedeffects::Vector{Int},randomeffects::Vector{Int},doseset::Tuple,nGstates::Int,nRsteps::Int,nalleles::Int,scheme::Int,totalmcmc::Int,total::Int,totalf::Int,tmax::Float64,nbursts::Int,nbchor::Int,temp=[1;1])

    # MCMC Metropolis-Hastings algorithm with truncated Normal proposal distribution to fit generalized telegraph model to live cell, mRNA FISH data, and burst correlations

    # n+1 is the number of Gstates, zeta is the number of R steps

    # rin is a (npr x totalsets) 2D array of the initial guess rate parameters

    # dataFISH is an array of FISH mRNA unnormalized count histograms
    # FISHsets specifies which columns of r pertain to FISH histograms

    # dataLC is an array of ON and OFF time unnormalized cummulative ON and OFF times
    # LCsets specifies which columns of r pertain to LC CDFs

    # fixed is the set of fixed effects params, i.e. complement of random effects
    # ignore is the set of params that are not fit
    # doseset is the set of data sets that have same dose and thus same random effects

    # piinv = 1/pi
    # pi2inv = 1/sqrt(pi)

    n = nGstates -1
    zeta = nRsteps

    println(n," ",zeta," ",randomeffects)

    nallelesf = nalleles
    ndoses = length(doseset)

    r = copy(rin)
    r = clamp.(r,1e-6,5000.)
    npr = length(r[:,1])
    totalsets = length(r[1,:])

    # set lower and upper bounds on initiation and decay rates
    upperbound = fill(5000.,npr)
    lowerbound = fill(1e-6,npr)
    lbvarF = 1e-3

    # total number of rate parameters
    reactiondecay = 2*n + zeta + 2
    if npr > reactiondecay
      rcorr = 1
    else
      rcorr = 0
    end

    nparams = 2*n + zeta + 2 + Int(rcorr>0)

    if nparams != npr
        print("Error, r has length ", npr,'\n')
        print("Should be ",nparams,'\n')
        return 0
    end

    # Make sure fixed effects are uniform across datasets
    # r[fixedeffects,:] = rin[fixedeffects,1]

    # Prepare FISH
    nFISH = length(dataFISH)

    if nFISH > 0
        # Declare FISH variables and prepare FISH histograms
        nhist = Array{Int}(nFISH)
        histFISH = Array{Array{Float64,1},1}(nFISH)
        countsFISH =Array{Int}(nFISH)

        for i = 1:nFISH
            nhist[i] = length(dataFISH[i][:,1])
            countsFISH[i] = sum(dataFISH[i][:,1])
            histFISH[i] = dataFISH[i][:,1]/countsFISH[i]
        end
    end

    # Prepare LC CDFs

    nLC = length(dataLC)
    if nLC > 0
        nbins, nCDF = size(dataLC[1][:,2:end])
        bins = dataLC[1][:,1]

        nCDF2 = div(nCDF,2)

        cdfON = Array{Array{Float64,1}}(nLC)
        cdfOFF = Array{Array{Float64,1}}(nLC)
        countsON = zeros(Int,nLC)
        countsOFF = zeros(Int,nLC)
        for kd = 1:nLC
            # kd = nDwell + k
            cdfON[kd] = dataLC[kd][:,2]
            cdfOFF[kd] = dataLC[kd][:,3]
            countsON[kd] = cdfON[kd][end]
            countsOFF[kd] = cdfOFF[kd][end]
            cdfON[kd] /= cdfON[kd][end]
            cdfOFF[kd] /= cdfOFF[kd][end]
        end
        Corrset = LCsets[1]
    else
        Corrset = FISHsets[3]
    end

    println(Corrset)

    # Check that r has correct length
    if totalsets != nFISH + nLC
        print("Error, r has ",totalsets," number of columns",'\n')
        print("Should be ",nFISH + nLC,'\n')
        return 0
    end

    # Compute initial Cost
    chiLC = 0
    chiFs = 0
    chiC = 0

    if nLC > 0
        # costLC(LCsets,cdfOFF,cdfON,countsOFF,countsON,bins,r,n,zet,total::Int,nbursts::Int,nallelesf::Int,scheme::Int,totalf::Int,tmax::Float64)
        chiLC = costLC(LCsets,cdfOFF,cdfON,countsOFF,countsON,bins,r,n,zeta,nallelesf,total,nbursts,scheme,totalf,tmax,temp[1])
    end
    if nFISH > 0
        chiFs = costFISH(FISHsets,histFISH,nhist,countsFISH,r,n,zeta,nallelesf,totalf,tmax,lbvarF,temp[2])
    end
    if rcorr > 0
        chiC = chicorr(nallelesf,n,zeta,r[:,Corrset],nbchor)
    end

    chi = chiLC + chiFs + chiC

    if isnan(chi)
        chi = 1e6
    end

    println("LC: ",chiLC," FISH: ",chiFs," Full: ",chi)

    # Set test parameters
    rt = copy(r)
    nr = length(r[:,1])

    # Set maximum likelihood variables
    rml = copy(r)
    chiml = chi
    chiFml = chiFs
    chiLCml = chiLC

    # Declare output variables
    chiout = Array{Float64}(totalmcmc+1)
    rout = Array{Array}(totalmcmc+1)
    chiout[1] = chi
    rout[1] = r

    pc = 0.05  # MCMC step size, fraction of current rate

    # MCMC loop
    for step = 1:totalmcmc

        if step == Int(div(totalmcmc,2))
            println("Halfway: ",now())
        end

        MHfactor = 1.

        # Select fixed effects parameters
        for i in fixedeffects
            d = Distributions.TruncatedNormal(r[i,1],pc*min(r[i,1],upperbound[i]),lowerbound[i],upperbound[i])
            rt[i,:] = repmat([rand(d)],1,totalsets)
            dt = Distributions.TruncatedNormal(rt[i,1],pc*min(rt[i,1],upperbound[i]),lowerbound[i],upperbound[i])
            # Compute MH nonsymmetric proposal distribution factor
            MHfactor *= Distributions.pdf(dt,r[i,1])/Distributions.pdf(d,rt[i,1])
        end

        # Select random effects parameters
        for k = 1:ndoses
            for j in randomeffects
                d = Distributions.TruncatedNormal(r[j,doseset[k][1]],pc*min(r[j,doseset[k][1]],1),lowerbound[j],upperbound[j])
                rt[j,doseset[k]] = rand(d)
            end
        end

        if nLC > 0
            chiLC = costLC(LCsets,cdfOFF,cdfON,countsOFF,countsON,bins,rt,n,zeta,nallelesf,total,nbursts,scheme,totalf,tmax,temp[1])
        end
        if nFISH > 0
            chiFs = costFISH(FISHsets,histFISH,nhist,countsFISH,rt,n,zeta,nallelesf,totalf,tmax,lbvarF,temp[2])
        end
        if rcorr > 0
            chiC = chicorr(nallelesf,n,zeta,rt[:,Corrset],nbchor)
        end

        chit = chiLC + chiFs + chiC

        if isnan(chit)
            chit = 1e6
        end

        if rand() < exp(.5*(chi - chit))*MHfactor
            chi = chit
            r = copy(rt)
            if chi < chiml
                chiml = chi
                chiLCml = chiLC
                chiFml = chiFs
                rml = copy(r)
            end
        end

        chiout[step+1] = chi
        rout[step+1] = r

    end

    rarray = zeros(totalmcmc+1,nparams,totalsets)

    for i = 1:totalmcmc+1
        rarray[i,:,:] = rout[i]
    end

    r025 = zeros(1,nparams,totalsets)
    r975 = similar(r025)
    for j = 1:nparams, k = 1:totalsets
        r025[1,j,k] = quantile(rarray[:,j,k],0.025)
        r975[1,j,k] = quantile(rarray[:,j,k],0.975)
    end

    rmean = mean(rarray,1)
    rstd = std(rarray,1)
    rmedian = median(rarray,1)

    return rml, [chiml mean(chiout) std(chiout) chiLCml chiFml], chiout, [rmean;rstd], [n zeta],[rmedian;r025;r975], rout[end]

end

function costLC(LCsets::Array{Int,1},cdfOFF::Array{Array{Float64,1},1},cdfON::Array{Array{Float64,1},1},countsOFF::Array{Int,1},countsON::Array{Int,1},bins::Array{Float64,1},rt::Array{Float64,2},n::Int,zet::Int,nallelesf::Int,total::Int,nbursts::Int,scheme::Int,totalf::Int,tmax::Float64,kT::Float64)
    # Live cell cost function
    chiLC = 0.
    for k in eachindex(cdfON)
        kr = LCsets[k]
        if scheme == 0
            modelOFF, modelON = @fastmath dwelltimeCDF(bins,n,zet,rt[:,kr],nallelesf,false)
        else
            modelOFF, modelON = @fastmath telegraphprefast(bins,n,zet,rt[:,kr],total,countsOFF[k],nallelesf)
        end
        modelOFF /= modelOFF[end]
        modelON /= modelON[end]
        DnOFF, indOFF = findmax(abs.(modelOFF-cdfOFF[k]))
        DnON, indON = findmax(abs.(modelON-cdfON[k]))
        chiLC += .5*countsOFF[k]*DnOFF^2 + .5*countsON[k]*DnON^2
    end
    chiLC *= 2/kT
end

function costFISH(FISHsets::Array{Int,1},histFISH::Array{Array{Float64,1},1},nhist::Array{Int64,1},countsFISH::Array{Int,1},rt::Array{Float64,2},n::Int,zet::Int,nallelesf::Int,totalf::Int,tmax::Float64,lbvarF::Float64,kT::Float64)
    # FISH cost function
    chiFs = 0.
    index0 = 1
    lbvar = 1e-5
    # lbarvar = 5./countsFISH[i]
    for i in eachindex(nhist)
        ir = FISHsets[i]
        histF = @fastmath telegraphprefast(n,zet,rt[:,ir],totalf,tmax,nhist[i],nallelesf,false)
        varF = kT*max.(histF.*(1-histF),lbvar)/countsFISH[i]*nhist[i]
        chiFs += sum((histFISH[i][index0:end]-histF[index0:end]).^2./varF[index0:end])
    end
    return chiFs
end
