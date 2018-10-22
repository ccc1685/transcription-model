# Function to compute ON/OFF time distributions for a single uncoupled allele

function dwelltimeSingleAlleleCDF(t::Vector,n::Int,na::Int,zeta::Int,rin::Vector)
	# Dwell time distribution of generalized telegraph model

	r = copy(rin)

	# decay reaction index
	reactiondecay = 2*n + zeta + 2

	nA = (n+1)*2^zeta
	nI = (n+1)*2

	# Declare Transition matrices
	T = zeros(nA,nA)
	TA = zeros(nA,nA)
	TI = zeros(nI,nI)
	SAinit = zeros(nA)
	SIinit = zeros(nI)

	mask = zeros(Int,n+1)

	ntime = length(t)

	# Declare survival and dwell time distributions
	SAj = zeros(ntime)
	SIj = zeros(ntime)

	SA = zeros(ntime,nA)
	SI = zeros(ntime,nI)

	gammap = zeros(n+2)
	gamman = zeros(n+2)
	nu = zeros(zeta+1)

	for i = 1:n
		gammap[i+1] = r[2*(i-1)+1]
		gamman[i+1] = r[2*i]
	end

	for i = 1:zeta+1
		nu[i] = r[2*n+i]
	end

	# Generate TA and T matrix
	for i=1:n+1, z=1:2^zeta, ip=1:n+1, w=1:2^zeta
		a = i + (n+1)*(z-1)
		b = ip + (n+1)*(w-1)
		zdigits = digits(z-1,2,zeta)
		wdigits = digits(w-1,2,zeta)
		z1 = zdigits[1]
		w1 = wdigits[1]
		zzeta = zdigits[zeta]
		wzeta = wdigits[zeta]
		zbar1 = zdigits[2:zeta]
		wbar1 = wdigits[2:zeta]
		zbarzeta = zdigits[1:zeta-1]
		wbarzeta = wdigits[1:zeta-1]
		if sum(zbarzeta)==0 && zzeta == 1
			mask[i] = a
		end
		TA[a,b] = (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w)*float(sum(wdigits)!=0) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)*(sum(wdigits)!=0)) + nu[zeta+1]*float((i==ip)*(zbarzeta==wbarzeta)*((zzeta==0)-(zzeta==1))*(wzeta==1))
		T[a,b] =  (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)) + nu[zeta+1]*float((i==ip)*(zbarzeta==wbarzeta)*((zzeta==0)-(zzeta==1))*(wzeta==1))
		for j = 1:zeta-1
			zbarj = zdigits[[1:j-1;j+2:zeta]]
			wbarj = wdigits[[1:j-1;j+2:zeta]]
			zj = zdigits[j]
			zj1 = zdigits[j+1]
			wj = wdigits[j]
			wj1 = wdigits[j+1]
			TA[a,b] += nu[j+1]*float((i==ip)*(zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
			T[a,b] += nu[j+1]*float((i==ip)*(zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
		end
		if zeta == 1
			SAinit[a] = float((z1==1)*(i==n+1))
		else
			SAinit[a] = float((z1==1)*(sum(zbar1)==0)*(i==n+1))
		end
	end

	# Generate TI matrix
	for i=1:n+1, z=0:1, ip=1:n+1, w=0:1
		a = i + (n+1)*z
		b = ip + (n+1)*w
		TI[a,b] = (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w)*float(w==0) + nu[1]*float((i==ip)*(i==n+1)*((z==1)-(z==0))*(w==0))
		SIinit[a] = float(z==0)
	end

	# Compute survival probability

	TAeiv = eigfact(TA);
	TAvects = TAeiv[:vectors]
	TAvals = TAeiv[:values]

	TAparams = TAvects\SAinit

	for j = 1:nA
		for i = 1:nA
			SAj += real(TAparams[i]*TAvects[j,i]*exp(TAvals[i]*t))
		end
		SA[:,j]=SAj
		SAj = zeros(ntime)
	end

	for i=1:n+1, z=1:2^zeta
		zdigits = digits(z-1,2,zeta)
		if sum(zdigits) == 0
			a = i + (n+1)*(z-1)
			SAj += SA[:,a]
		end
	end

	TIeiv = eigfact(TI);
	TIvects = TIeiv[:vectors]
	TIvals = TIeiv[:values]

	pss = nullspace(T)
	y = [pss[mask]/sum(pss[mask]);zeros(div(length(SIinit),2))]
	SIinit = SIinit.*y
	TIparams = TIvects\SIinit

	for j = 1:nI
		for i = 1:nI
			SIj += real(TIparams[i]*TIvects[j,i]*exp(TIvals[i]*t))
		end
		SI[:,j]=SIj
		SIj = zeros(ntime)
	end

	SIj = SI[:,nI]
	SIj /= SIj[end]
	SAj /= SAj[end]

	#return DA,DAj,SA,SAj,SAinit,DI,DIj,SI,SIj,SIinit,TA,TI,T,mask
	return SIj,SAj

end
