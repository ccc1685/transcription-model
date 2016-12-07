


function speedtest()
# Dwell time PDF of generalized telegraph model

n = 2
zeta = 3
rang = 4000
nalleles = 5
#total = 1000000000
total = 1000000
dt = .1
rin = [5.38397e-5  0.00128317  0.0277282  2.86347  1.23516  0.205716  7.22633  0.146527  0.00076717  1.36376]

# copy rates
r = copy(rin)

# decay reaction index
reactiondecay = 2*n + zeta + 2

# histogram array
nhist = Int(cld(rang,dt))+1
#histofftdd = zeros(Int,nhist)

# time
t=0.

# tIA and TAI are times when gene state transitions between active and inactive states
# they are used to compute dwell times
tIA = zeros(Float64,nalleles)
tAI = zeros(Float64,nalleles)

# active start site
mu = n
#print(n," ",mu," ",zeta,'\n')

# Initialize rates
gammaf = zeros(n+1)
gammab = zeros(n+1)
nu = Array(Float64,zeta+1)

# assign gene rates
for i = 1:n
    gammaf[i] = r[2*(i-1)+1]
    gammab[i+1] = r[2*i]
end

# assign pre-mRNA rates
for i = 1:zeta+1
  nu[i] = r[2*n + i]
end

# assign decay rate
delta = r[reactiondecay]

# assign correlation rate
if length(r) > reactiondecay
  rcorr = r[reactiondecay+1]
else
  rcorr = 0
end

# Initialize to gene state 1 to 1, all other states are zero
state = ones(Int,nalleles)
nindex = n + 1

# pre-mRNA states
z = Array(Array,nalleles)

# initialize propensities
af = Array(Float64,nalleles)
ab = similar(af)
apre1 = similar(af)
az = Array(Array,nalleles)

apremid = zeros(nalleles)     # sum(az[2:zeta]) middle pre-mRNAstep
afree = zeros(nalleles)       # az[zeta+1]  create free mRNA

for i = 1:nalleles
    z[i] = zeros(Int8,zeta)  # Intialize pre-mRNA state to all 0
    af[i] = gammaf[state[i]]
    ab[i] = gammab[state[i]]
    az[i] = zeros(zeta+1)
    apre1[i] = float(1-z[i][1])*nu[1]*float(state[i]==n+1)   #first pre-mRNA step
end
astate = af + ab

#adecay = 0       #m*delta

amrna = apre1 + apremid + afree

for steps = 1:total

asum = sum(amrna) + sum(astate) #+ adecay

if asum > 0
# choose which state to update and time of update
   r = asum*rand()
   tau = -log(rand())/asum

# update time
   t += tau
   #print(tau,' ',asum,' ',state,'\n')

   #times[steps]=tau
   #mRNA[steps] = m
#=
   if r > sum(amrna + astate)
   # mRNA decay
         m -= 1
      	 adecay = m*delta

   else
=#
    agenecs = cumsum(amrna + astate)
      l = 1
      while r > agenecs[l]
    	 l += 1
      end
      if l > 1
           r -= agenecs[l-1]
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

      for k in collect(1:nalleles)
	        if state[k] == mu && state[l] == mu + 1
	          af[k] = gammaf[state[k]] + rcorr
	          astate[k] = af[k] + ab[k]
          else
            af[k] = gammaf[state[k]]
	          ab[k] = gammab[state[k]]
	          astate[k] = af[k] + ab[k]
	      end
      end

      # update propensities
      af[l] = gammaf[state[l]]
      ab[l] = gammab[state[l]]
      astate[l] = af[l] + ab[l]

      if state[l] > mu
       	 apre1[l] = float(1-z[l][1])*nu[1]
       	# apre1 = az[l][1]
      else
         apre1[l] = 0
      end

      amrna[l] = apremid[l] + apre1[l] + afree[l]

     else

     # premRNA update
    if r > apremid[l] + apre1[l]
        # increase free mRNA
        #m += 1
	      #adecay = m*delta

	      z[l][zeta] = 0
	          #az[l][zeta+1] = 0
	      afree[l] = 0

	      if sum(z[l]) == 0
           tAI[l] = t
#           tactive[acount] = tAI[l] - tIA[l];
#           acount +=1
	        end

	      if zeta == 1
	    if state[l] > mu
	       #az[l][1] = nu[1]
	       apre1[l] = nu[1]
	    else
		      #az[l][1] = 0
		        apre1[l] = 0
	    end
	 else
	    az[l][zeta] = float(z[l][zeta-1])*nu[zeta]
	    apremid[l] = sum(az[l][2:zeta])
	 end

	 amrna[l] = apre1[l] + apremid[l] + afree[l]
   elseif r > apremid[l]
       # increase first pre-mRNA state

         if sum(z[l]) == 0
           tIA[l] = t
	         tinactive = ceil(Int,(tIA[l] - tAI[l])/dt)
           if tinactive <= nhist && tinactive > 0
             #histofftdd[tinactive] += 1
             #histofftdd[1] = 1
           end
         end

         z[l][1] = 1
         #az[l][1] = 0
	 apre1[l] = 0

	 if zeta == 1
	    #az[l][2] = nu[2]
	    afree[l] = nu[2]
	 else
	    az[l][2] = float(1-z[l][2])*nu[2]
	    apremid[l] = sum(az[l][2:zeta])
	 end

         amrna[l] = apre1[l] + apremid[l] + afree[l]

       else
       # update mid pre-mRNA state

         azsum = cumsum(az[l][2:zeta])
	 k = 1
	 while r > azsum[k]
	       k += 1
	 end
	 i = k + 1

	 z[l][i] = 1
         z[l][i-1] = 0
	 az[l][i] = 0

	 if i == zeta
	    #az[l][i+1] = nu[i+1]
	    afree[l] = nu[i+1]
	 else
	    az[l][i+1] = float(1-z[l][i+1])*nu[i+1]
	 end
	 if i == 2
	   if state[l] > mu
	      #az[l][1] = nu[1]
	      apre1[l] = nu[1]
	   else
	      az[l][1] = 0
              apre1[l] = 0
	   end
	 else
	    az[l][i-1] = float(z[l][i-2])*nu[i-1]
	 end
	 apremid[l] = sum(az[l][2:zeta])
	 amrna[l] = apre1[l] + apremid[l] + afree[l]

       end # if r > apremid[l] + apre1[l]
    end # if r > amrna[l]
   #end #if r > sum(amrna + adecay)

end  #if asum > 0

end  # for steps = 1:total

#return histofftdd

end
