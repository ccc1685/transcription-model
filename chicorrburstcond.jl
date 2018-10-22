# Functions to compute burst-to-burst allele cross_correlations and burst number correlations


function chicorr(nalleles::Int,n::Int,zet::Int,rin::Vector,nbursts::Int)
  # Number of burst simulations for nalleles alleles in same cell with correlation between initial conditions of bcor

  # Bounds to prevent infinite loop
  nexps = 10000
  nsteps = 50000

  start0 = ones(Int,nalleles)

  bursts = zeros(Int,nexps,nalleles)
  gfprow = Array{Int,1}(nalleles)
  times = Array{Float64,1}(nsteps)
  gfp = Array{Int,2}(nsteps,nalleles)

  total = 854 # ~ 14 hours

  # start0[1] = 3
  start0[2] = 2

  # active start site
  mu = n

  n,zet,nalleles,gammaf,gammab,nu,rcorr,state,z,af,ab,az,apre1,apremid,afree = telegraphcorr(n,zet,rin,start0,nalleles)

  astate = af + ab
  amrna = apre1 + apremid + afree

  cc = 0.
  avg = 0.
  totaltime = 0.
  iexps = 1
  ibursts = 0

  while ibursts < nbursts && iexps <= nexps

    t = 0
    steps = 0

    while t <= total && steps < nsteps

      asum = sum(amrna + astate)

      steps += 1

      # update gfp
      for k = 1:nalleles
        gfprow[k] = sum(z[k])
      end
      gfp[steps,:] = gfprow

      # choose which state to update and time of update
      r = asum*rand()
      tau = -log(rand())/asum

      # update time
      t += tau
      times[steps] = t

      agenecs = amrna[1] + astate[1]
      agene0 = 0.
      l = 1
      while r > agenecs
        agene0 = agenecs
        agenecs += amrna[1+l] + astate[l+1]
        l += 1
      end
      if l > 1
        r -= agene0
      end

      if r > amrna[l]

        # state update
        if r > ab[l] + amrna[l]
          # update forward reaction
          state[l] += 1
        else
          # update backward reaction
          state[l] -= 1
        end

        # update propensities
        af[l] = gammaf[state[l]]
        ab[l] = gammab[state[l]]
        astate[l] = af[l] + ab[l]

        # Activate coupling if any alleles are in state mu + 1
        ba = false
        for k = 1:nalleles
          ba |= state[k] > mu
        end
        if ba
          for k = 1:nalleles
            if state[k] == mu
              af[k] = gammaf[state[k]] + rcorr
              astate[k] = af[k] + ab[k]
            end
          end
        else
          for k = 1:nalleles
            if state[k] == mu
              af[k] = gammaf[state[k]]
              astate[k] = af[k] + ab[k]
            end
          end
        end

        if state[l] > mu
          apre1[l] = float(1-z[l][1])*nu[1]
        else
          apre1[l] = 0.
        end

        amrna[l] = apremid[l] + apre1[l] + afree[l]

      else

        # premRNA update
        if r > apremid[l] + apre1[l]

          z[l][zet] = 0
          afree[l] = 0.

          if zet == 1
            if state[l] > mu
              apre1[l] = nu[1]
            else
              apre1[l] = 0.
            end
          else
            az[l][zet] = float(z[l][zet-1])*nu[zet]
            apremid[l] = sum(az[l][2:zet])
          end

          amrna[l] = apre1[l] + apremid[l] + afree[l]

        elseif r > apremid[l]
          # increase first pre-mRNA state

          if sum(z[l]) == 0
            bursts[iexps,l] += 1
          end

          z[l][1] = 1
          apre1[l] = 0.

          if zet == 1
            afree[l] = nu[2]
          else
            az[l][2] = float(1-z[l][2])*nu[2]
            apremid[l] = sum(az[l][2:zet])
          end

          amrna[l] = apre1[l] + apremid[l] + afree[l]

        else
          # update mid pre-mRNA state
          azsum = az[l][2]
          k = 1
          while r > azsum
            azsum += az[l][2+k]
            k += 1
          end

          i = k + 1

          z[l][i] = 1
          z[l][i-1] = 0
          az[l][i] = 0.

          if i == zet
            afree[l] = nu[i+1]
          else
            az[l][i+1] = float(1-z[l][i+1])*nu[i+1]
          end
          if i == 2
            if state[l] > mu
              apre1[l] = nu[1]
            else
              az[l][1] = 0.
              apre1[l] = 0.
            end
          else
            az[l][i-1] = float(z[l][i-2])*nu[i-1]
          end
          apremid[l] = sum(az[l][2:zet])
          amrna[l] = apre1[l] + apremid[l] + afree[l]

        end # if r > apremid[l] + apre1[l]
      end # if r > amrna[l]
    end  # while t <= total

    times[steps] = total
    # println(sum(bursts[i,:]))

    nt = length(times[1:steps])
    G = zeros(nt,nalleles)
    #mins = round(Int,times[i])+1
    nactive = 0
    k = 0
    for j = 1:nalleles
      if maximum(gfp[1:steps,j]) > .5
        nactive += 1
        k += 1
        G[:,k] = gfp[1:steps,j]
      end
    end
    if nactive > 1
      dt = diff(times[1:steps])
      dt = [times[1];dt]
      for j = 1:nactive-1, k = j+1:nactive
        cc += sum(G[:,j].*G[:,k].*dt)
      end
      totaltime += nactive*(nactive-1)
      avg += sum(G.*dt)*(nactive-1)
    end
    # ibursts = sum(sum(bursts.>0,2).>1)
    ibursts += length(find(bursts[iexps,:])) > 1
    iexps += 1
  end #while ibursts < nbursts && iexps < nexps
  # println("corr:",ibursts)
  # println(iexps)

  totaltime *= total

  avg /= totaltime

  cc /= totaltime
  cc *= 2
  cc /= avg^2
  cc -= 1

  if isnan(cc)
      cc = -10.
  end

  cb = Array{Float64}(nexps)

  nalleles = size(bursts)[2]

  # burst2 = Array{Int,2}(0,2)
  # for i = 1:nalleles-1, j = i+1:nalleles
  #   burst2 = vcat(burst2,[bursts[:,i] bursts[:,j]])
  # end
  npairs = div(nalleles*(nalleles-1),2)
  burst2 = zeros(Int,nexps*npairs,2)
  count = 1
  for i = 1:nalleles-1, j = i+1:nalleles
    burst2[(count-1)*nexps+1:count*nexps,1] = bursts[:,i]
    burst2[(count-1)*nexps+1:count*nexps,2] = bursts[:,j]
    count += 1
  end

  burst2 = burst2[burst2[:,1].*burst2[:,2].>0,:]
  if isempty(burst2[:,1]) || isempty(burst2[:,1])
      cb = -10.
  else
      cb = cor(burst2[:,1],burst2[:,2])
  end

  if isnan(cb)
      cb = -10.
  end

  # cbo = bootstrapburst(bursts,76,2*nexps)
  # cb = mean(cbo)

  # println(cb," ",std(cbo))
  # println(cc, " ",cb)

  sigb = .1
  sigcc = .1875

  # sigb = .2
  # sigcc = .2

  # return exp(-ac*(cb -.486)^2-am*(cc-0.5)^2)
  return (cb -.486)^2/(sigb^2) + (cc-0.4202)^2/(sigcc^2)

end # function

function telegraphcorr(n::Int,zet::Int,r::Vector,state::Vector,nalleles::Int)

	# Generalized n,mu,zet, telegraph model,
	# n+1 total gene states
	# gene states are labeled from 1 to n+1
	# zet pre-mRNA steps
	# There are 2n gene reactions, na*zet pre-mRNA reactions, and 1 mRNA decay reaction
	# rin = reaction rates
	# forward gene reactions are labeled by odd numbers <2n, backward reactions are even numbers <= 2n
	# reactions [2n+1, 2n+zet] are pre-mRNA recations,organized by starting spot, reaction 2n+(zet+1)+1 is mRNA decay reaction

	# Use Gillespie direct method with speed up from Gibson and Bruck

	# decay reaction index
	reactiondecay = 2*n + zet + 2

	# Initialize rates
	gammaf = zeros(n+1)
	gammab = zeros(n+1)
	nu = Array{Float64,1}(zet+1)

	# assign gene rates
	for i = 1:n
		gammaf[i] = r[2*(i-1)+1]
		gammab[i+1] = r[2*i]
	end

	# assign pre-mRNA rates
	for i = 1:zet+1
		nu[i] = r[2*n + i]
	end

	# assign decay rate
	delta = r[reactiondecay]

	# assign correlation rate
	if length(r) > reactiondecay
		rcorr = r[reactiondecay+1]
	else
		rcorr = 0.
	end

	# Initialize gene state to state0
	# state = copy(state0)

	# pre-mRNA states
	z = Array{Array{Int,1},1}(nalleles)

	# initialize propensities
	af = Array{Float64,1}(nalleles)
	ab = Array{Float64,1}(nalleles)
	apre1 = Array{Float64,1}(nalleles)
	az = Array{Array{Float64,1}}(nalleles)

	apremid = zeros(nalleles)     # sum(az[2:zet]) middle pre-mRNAstep
	afree = zeros(nalleles)       # az[zet+1]  create free mRNA

	for i = 1:nalleles
		z[i] = zeros(Int8,zet)  # Intialize pre-mRNA state to all 0
		af[i] = gammaf[state[i]]
		ab[i] = gammab[state[i]]
		az[i] = zeros(zet+1)
		apre1[i] = float(1-z[i][1])*nu[1]*float(state[i]==n+1)   #first pre-mRNA step
	end

	astate = af + ab

	amrna = apre1 + apremid + afree

	return n,zet,nalleles,gammaf,gammab,nu,rcorr,state,z,af,ab,az,apre1,apremid,afree

end

function chicorr(nalleles::Int,n::Int,zet::Int,rin::Vector)

  # nbursts = 760
  nbursts = 500

  chicorr(nalleles,n,zet,rin,nbursts)
end
