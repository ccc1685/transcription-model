#Functions for computing ON/OFF LC histograms using first passage time solution of Master equation


function dwelltimeCDF(t::Vector,n::Int,na::Int,zeta::Int,rin::Vector)
	dwelltimeCDF(t,n,zeta,rin,5,false)
end


function dwelltimeCDF(t::Vector,n::Int,zeta::Int,rin::Vector,nallelesin::Int,ss::Bool)
	# Dwell time distribution of a single allele in a multi-allele model

	r = copy(rin)

	# decay reaction index
	reactiondecay = 2*n + zeta + 2

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

	# assign correlation rate
	if length(r) > reactiondecay
		rcorr = r[reactiondecay+1]
		nalleles = nallelesin
	else
		rcorr = 0.
		nalleles = 1
	end

	# println(rcorr)
	nG = n+1
	nA = (nG)*2^zeta
	nI = (nG)*2


	# Declare Transition matrices
	T = zeros(nA,nA)
	TA = zeros(nA,nA)
	TI = zeros(nI,nI)
	TG = zeros(nG,nG)

	# Initial conditions
	SAinit = zeros(nA)
	SIinit = zeros(nI*nG^(nalleles-1))

	# filters
	mask = zeros(Int,nG)
	maskA = zeros(Int,nA)
	maskzb0z1 = fill(false,nA)
	maskIinit = fill(false,nI)

	# Coupling matrices
	XYA = fill(false,nA,nA)
	YYA = fill(false,nA,nA)
	XXA = fill(false,nA,nA)

	# Generate TA and T transition matrices for a single allele
	for i=1:nG, z=1:2^zeta
		for ip=1:nG, w=1:2^zeta
			a = i + nG*(z-1)
			b = ip + nG*(w-1)
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
				# Needed for SIinit
				mask[i] = a
				maskzb0z1[a] = true
			end
			if sum(zdigits) == 0
				maskA[a] = 1
			end
			if i == n && ip == n && z == w
				YYA[a,b] = true
			end
			if i == nG && ip == n && z == w
				XYA[a,b] = true
			end
			if i == nG && ip == nG && z == w
				XXA[a,b] = true
			end
			TA[a,b] = (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w)*float(sum(wdigits)!=0) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)*(sum(wdigits)!=0)) + nu[zeta+1]*float((i==ip)*(zbarzeta==wbarzeta)*((zzeta==0)-(zzeta==1))*(wzeta==1))
			T[a,b] =  (gammap[i]*(i-1==ip) - (gamman[i]+gammap[i+1])*(i==ip) + gamman[i+1]*(i+1==ip))*(z==w) + nu[1]*((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)) + nu[zeta+1]*((i==ip)*(zbarzeta==wbarzeta)*((zzeta==0)-(zzeta==1))*(wzeta==1))
			for j = 1:zeta-1
				zbarj = zdigits[[1:j-1;j+2:zeta]]
				wbarj = wdigits[[1:j-1;j+2:zeta]]
				zj = zdigits[j]
				zj1 = zdigits[j+1]
				wj = wdigits[j]
				wj1 = wdigits[j+1]
				TA[a,b] += nu[j+1]*float((i==ip)*(zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
				T[a,b] += nu[j+1]*((i==ip)*(zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
			end
			if zeta == 1
				SAinit[a] = float((z1==1)*(i==n+1))
			else
				SAinit[a] = float((z1==1)*(sum(zbar1)==0)*(i==n+1))
			end
		end
	end

	XYI = fill(false,nI,nI)
	YYI = fill(false,nI,nI)
	XXI = fill(false,nI,nI)

	# Generate TI transition matrix for a single allele
	for i=1:n+1, z=0:1
		for ip=1:n+1, w=0:1
			a = i + (n+1)*z
			b = ip + (n+1)*w
			if i == n && ip == n && z == w
				YYI[a,b] = true
			end
			if i == n+1 && ip == n  && z == w
				XYI[a,b] = true
			end
			if i == n+1 && ip == n+1 && z == w
				XXI[a,b] = true
			end
			TI[a,b] = (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w)*float(w==0) + nu[1]*float((i==ip)*(i==n+1)*((z==1)-(z==0))*(w==0))
			maskIinit[a] = float(z==0)
		end
	end

	XYG = fill(false,nG,nG)
	YYG = fill(false,nG,nG)
	XXG = fill(false,nG,nG)

	# Generate transition matrices for G states of a single allele
	for i=1:nG
		for ip=1:nG
			a = i
			b = ip
			if i == n && ip == n
				YYG[a,b] = true
			end
			if i == nG && ip == n
				XYG[a,b] = true
			end
			if i == nG && ip == nG
				XXG[a,b] = true
			end
			TG[a,b] =  gammap[i]*(i-1==ip) - (gamman[i]+gammap[i+1])*(i==ip) + gamman[i+1]*(i+1==ip)
		end
	end

	# nT = nA^nalleles

	VA = Array{Float64,2}
	VI = Array{Float64,2}
	VG = Array{Float64,2}
	V = Array{Float64,2}

	maskI = zeros(Int,nI)
	maskI[end] = 1

	if nalleles > 1
		# Assemble full transition matrix

		I = Array{Array{Int64,2},1}(nalleles-1)
		II = Array{Array{Int64,2},1}(nalleles-1)
		IA = Array{Array{Int64,2},1}(nalleles-1)

		I[1] = eye(Int,nG)
		II[1] = eye(Int,nI)
		IA[1] = eye(Int,nA)
		for i = 2:nalleles-1
			I[i] = eye(Int,nG^i)
			II[i] = eye(Int,nI*nG^(i-1))
			IA[i] = eye(Int,nA*nG^(i-1))
		end

		VI =  kron(II[nalleles-1],TG)
		VA =  kron(IA[nalleles-1],TG)

		# Assemble uncoupled full transition matrix using Kronecker product
		for j = 2:nalleles-1
			VI += kron(II[j-1],kron(TG,I[nalleles-j]))
			VA += kron(IA[j-1],kron(TG,I[nalleles-j]))
		end

		V = kron(T,I[nalleles-1]) + VA
		VI += kron(TI,I[nalleles-1])
		VA += kron(TA,I[nalleles-1])

		# Assign coupling

		if nalleles > 2
			Beye = Array{Array{Bool,2},1}(nalleles-2)
			BeyeI = Array{Array{Bool,2},1}(nalleles-2)
			BeyeA = Array{Array{Bool,2},1}(nalleles-2)
			for i = 1:nalleles-2
				Beye[i] = eye(Bool,nG^i)
				BeyeI[i] = eye(Bool,nI*nG^(i-1))
				BeyeA[i] = eye(Bool,nA*nG^(i-1))
			end
		end

		XB = Array{Array{Bool,2},1}(nalleles-1)
		XBI = Array{Array{Bool,2},1}(nalleles-1)
		XBA = Array{Array{Bool,2},1}(nalleles-1)

		XBA[1] = XXA
		XBI[1] = XXI
		XB[1] = XXG
		for j = 2:nalleles-1
			XB[j] = kron(Beye[j-1],XB[1]) .| kron(XB[1],Beye[j-1])
			XBA[j] = kron(BeyeA[j-1],XB[1]) .| kron(XBA[1],Beye[j-1])
			XBI[j] = kron(BeyeI[j-1],XB[1]) .| kron(XBI[1],Beye[j-1])
			for k = 1:j-2
				XB[j] .|= kron(kron(Beye[k],XXG),Beye[j-1-k])
				XBA[j] .|= kron(kron(BeyeA[k],XXG),Beye[j-1-k])
				XBI[j] .|= kron(kron(BeyeI[k],XXG),Beye[j-1-k])
			end
		end

		XYVA = kron(XYA,XB[nalleles-1])
		YYVA = kron(YYA,XB[nalleles-1])
		XYVI = kron(XYI,XB[nalleles-1])
		YYVI = kron(YYI,XB[nalleles-1])

		VI[XYVI] += rcorr
		VI[YYVI] -= rcorr

		VA[XYVA] += rcorr
		VA[YYVA] -= rcorr

		V[XYVA] += rcorr
		V[YYVA] -= rcorr

		for i = 1:nalleles-2

			XYVA = kron(BeyeA[i],kron(XYG,XB[nalleles-1-i])) .| kron(XBA[i],kron(XYG,Beye[nalleles-1-i]))
			YYVA = kron(BeyeA[i],kron(YYG,XB[nalleles-1-i])) .| kron(XBA[i],kron(YYG,Beye[nalleles-1-i]))
			XYVI = kron(BeyeI[i],kron(XYG,XB[nalleles-1-i])) .| kron(XBI[i],kron(XYG,Beye[nalleles-1-i]))
			YYVI = kron(BeyeI[i],kron(YYG,XB[nalleles-1-i])) .| kron(XBI[i],kron(YYG,Beye[nalleles-1-i]))

			VI[XYVI] += rcorr
			VI[YYVI] -= rcorr

			VA[XYVA] += rcorr
			VA[YYVA] -= rcorr

			V[XYVA] += rcorr
			V[YYVA] -= rcorr
		end

		XYVA = kron(XBA[nalleles-1],XYG)
		YYVA = kron(XBA[nalleles-1],YYG)
		XYVI = kron(XBI[nalleles-1],XYG)
		YYVI = kron(XBI[nalleles-1],YYG)

		VI[XYVI] += rcorr
		VI[YYVI] -= rcorr

		VA[XYVA] += rcorr
		VA[YYVA] -= rcorr

		V[XYVA] += rcorr
		V[YYVA] -= rcorr

		# Compute initial condition for Inactive state distribution
		# Steady state distribution is null space of V
		pss = nullspace(V)
		pss /= sum(pss)

		# Marginalize to obtain steady state distribution of 1 allele and nalleles-1 alleles
		nAnm1 = nG^(nalleles-1)
		pssnm1 = Array{Float64,1}(nAnm1)


		for i = 1:nAnm1
			pssnm1[i] = sum(pss[(i-1)*nA+1:i*nA])
		end

		if ss
			pss1 = Array{Float64,1}(nA)
			for i = 1:nA
				pss1[i] = sum(pss[(i-1)*nAnm1+1:i*nAnm1])
			end
			Gss = zeros(n+1)
			for i=1:n+1, z=1:2^zeta
				a = i + (n+1)*(z-1)
				Gss[i] += pss1[a]
			end
			println(Gss)
			return Gss
		end

		y = pss[kron(maskzb0z1,fill(true,nG^(nalleles-1)))]

		SIinit[kron(maskIinit,fill(true,nG^(nalleles-1)))] = y
		SIinit /= sum(SIinit)

		maskI = kron(maskI,ones(Int,nG^(nalleles-1)))

		SAinit = kron(SAinit,pssnm1)
		maskA = kron(maskA,ones(Int,nG^(nalleles-1)))


	else
		VA = TA
		VI = TI

		pss = nullspace(T)
		pss = pss/sum(pss)

		if ss
			Gss = zeros(n+1)
			for i=1:n+1, z=1:2^zeta
				a = i + (n+1)*(z-1)
				Gss[i] += pss[a]
			end

			println(Gss)
			return Gss
		end

		SIinit = [pss[mask]/sum(pss[mask]);zeros(div(length(SIinit),2))]

	end

	# Compute survival probability

	# Declare survival and dwell time distributions
	ntime = length(t)
	SAj = zeros(ntime)
	SIj = zeros(ntime)

	TAeiv = eigfact(VA);
	TAvects = TAeiv[:vectors]
	TAvals = TAeiv[:values]

	TAparams = TAvects\SAinit

	for j in find(maskA)
		for i = 1:size(TAvects,2)
			SAj += real(TAparams[i]*TAvects[j,i]*exp.(TAvals[i]*t))
		end
	end

	TIeiv = eigfact(VI)
	TIvects = TIeiv[:vectors]
	TIvals = TIeiv[:values]

	TIparams = TIvects\SIinit

	for j in find(maskI)
		for i = 1:size(TIvects,2)
			SIj += real(TIparams[i]*TIvects[j,i]*exp.(TIvals[i]*t))
		end
	end

	SIj /= SIj[end]
	SAj /= SAj[end]

	return SIj,SAj

end

function dwelltimeCoupledOFF(t::Vector,n::Int,zeta::Int,rin::Vector,nalleles::Int)
	# OFF time distribution of a single allele in a multi-allele model

	r = copy(rin)

	# decay reaction index
	reactiondecay = 2*n + zeta + 2

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

	# assign correlation rate
	if length(r) > reactiondecay
		rcorr = r[reactiondecay+1]
	else
		rcorr = 0.
	end

	# println(rcorr)
	nG = n+1
	nA = (nG)*2^zeta
	nI = (nG)*2

	# nT = nA

	# Declare Transition matrices
	T = zeros(nA,nA)
	TA = zeros(nA,nA)
	TI = zeros(nI,nI)
	TG = zeros(nG,nG)

	# Initial conditions
	SAinit = zeros(nA)
	SIinit = zeros(nI*nG^(nalleles-1))

	# filters
	mask = zeros(Int,nG)
	maskA = zeros(Int,nA)
	maskzb0z1 = fill(false,nA)
	maskIinit = fill(false,nI)

	# Coupling matrices
	XYA = fill(false,nA,nA)
	YYA = fill(false,nA,nA)
	XXA = fill(false,nA,nA)

	# Generate TA and T transition matrices for a single allele
	for i=1:nG, z=1:2^zeta
		for ip=1:nG, w=1:2^zeta
			a = i + nG*(z-1)
			b = ip + nG*(w-1)
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
				# Needed for SIinit
				mask[i] = a
				maskzb0z1[a] = true
			end
			if sum(zdigits) == 0
				maskA[a] = 1
			end
			if i == n && ip == n && z == w
				YYA[a,b] = true
			end
			if i == nG && ip == n && z == w
				XYA[a,b] = true
			end
			if i == nG && ip == nG && z == w
				XXA[a,b] = true
			end
			TA[a,b] = (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w)*float(sum(wdigits)!=0) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)*(sum(wdigits)!=0)) + nu[zeta+1]*float((i==ip)*(zbarzeta==wbarzeta)*((zzeta==0)-(zzeta==1))*(wzeta==1))
			T[a,b] =  (gammap[i]*(i-1==ip) - (gamman[i]+gammap[i+1])*(i==ip) + gamman[i+1]*(i+1==ip))*(z==w) + nu[1]*((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)) + nu[zeta+1]*((i==ip)*(zbarzeta==wbarzeta)*((zzeta==0)-(zzeta==1))*(wzeta==1))
			for j = 1:zeta-1
				zbarj = zdigits[[1:j-1;j+2:zeta]]
				wbarj = wdigits[[1:j-1;j+2:zeta]]
				zj = zdigits[j]
				zj1 = zdigits[j+1]
				wj = wdigits[j]
				wj1 = wdigits[j+1]
				TA[a,b] += nu[j+1]*float((i==ip)*(zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
				T[a,b] += nu[j+1]*((i==ip)*(zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
			end
			if zeta == 1
				SAinit[a] = float((z1==1)*(i==n+1))
			else
				SAinit[a] = float((z1==1)*(sum(zbar1)==0)*(i==n+1))
			end
		end
	end

	XYI = fill(false,nI,nI)
	YYI = fill(false,nI,nI)
	XXI = fill(false,nI,nI)
	# Generate TI transition matrix for a single allele
	for i=1:n+1, z=0:1
		for ip=1:n+1, w=0:1
			a = i + (n+1)*z
			b = ip + (n+1)*w
			if i == n && ip == n && z == w
				YYI[a,b] = true
			end
			if i == n+1 && ip == n  && z == w
				XYI[a,b] = true
			end
			if i == n+1 && ip == n+1 && z == w
				XXI[a,b] = true
			end
			TI[a,b] = (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w)*float(w==0) + nu[1]*float((i==ip)*(i==n+1)*((z==1)-(z==0))*(w==0))
			maskIinit[a] = float(z==0)
		end
	end

	XYG = fill(false,nG,nG)
	YYG = fill(false,nG,nG)
	XXG = fill(false,nG,nG)

	# Generate transition matrices for G states of a single allele
	for i=1:nG
		for ip=1:nG
			a = i
			b = ip
			if i == n && ip == n
				YYG[a,b] = true
			end
			if i == nG && ip == n
				XYG[a,b] = true
			end
			if i == nG && ip == nG
				XXG[a,b] = true
			end
			TG[a,b] =  gammap[i]*(i-1==ip) - (gamman[i]+gammap[i+1])*(i==ip) + gamman[i+1]*(i+1==ip)
		end
	end

	# nT = nA^nalleles

	VA = Array{Float64,2}
	VI = Array{Float64,2}
	VG = Array{Float64,2}
	V = Array{Float64,2}

	maskI = zeros(Int,nI)
	maskI[end] = 1

	if nalleles > 1

		I = Array{Array{Int64,2},1}(nalleles-1)
		II = Array{Array{Int64,2},1}(nalleles-1)
		IA = Array{Array{Int64,2},1}(nalleles-1)

		I[1] = eye(Int,nG)
		II[1] = eye(Int,nI)
		IA[1] = eye(Int,nA)
		for i = 2:nalleles-1
			I[i] = eye(Int,nG^i)
			II[i] = eye(Int,nI*nG^(i-1))
			IA[i] = eye(Int,nA*nG^(i-1))
		end

		VI =  kron(II[nalleles-1],TG)
		VA =  kron(IA[nalleles-1],TG)

		# Assemble uncoupled full transition matrix using Kronecker product
		for j = 2:nalleles-1
			VI += kron(II[j-1],kron(TG,I[nalleles-j]))
			VA += kron(IA[j-1],kron(TG,I[nalleles-j]))
		end

		V = kron(T,I[nalleles-1]) + VA
		VI += kron(TI,I[nalleles-1])
		VA += kron(TA,I[nalleles-1])

		# println(size(VA))
		# println(size(VI))

		# Assign coupling

		if nalleles > 2
			Beye = Array{Array{Bool,2},1}(nalleles-2)
			BeyeI = Array{Array{Bool,2},1}(nalleles-2)
			BeyeA = Array{Array{Bool,2},1}(nalleles-2)
			for i = 1:nalleles-2
				Beye[i] = eye(Bool,nG^i)
				BeyeI[i] = eye(Bool,nI*nG^(i-1))
				BeyeA[i] = eye(Bool,nA*nG^(i-1))
			end
		end

		XB = Array{Array{Bool,2},1}(nalleles-1)
		XBI = Array{Array{Bool,2},1}(nalleles-1)
		XBA = Array{Array{Bool,2},1}(nalleles-1)

		XBA[1] = XXA
		XBI[1] = XXI
		XB[1] = XXG
		for j = 2:nalleles-1
			XB[j] = kron(Beye[j-1],XB[1]) | kron(XB[1],Beye[j-1])
			XBA[j] = kron(BeyeA[j-1],XB[1]) | kron(XBA[1],Beye[j-1])
			XBI[j] = kron(BeyeI[j-1],XB[1]) | kron(XBI[1],Beye[j-1])
			for k = 1:j-2
				XB[j] |= kron(kron(Beye[k],XXG),Beye[j-1-k])
				XBA[j] |= kron(kron(BeyeA[k],XXG),Beye[j-1-k])
				XBI[j] |= kron(kron(BeyeI[k],XXG),Beye[j-1-k])
			end
		end

		XYVA = kron(XYA,XB[nalleles-1])
		YYVA = kron(YYA,XB[nalleles-1])
		XYVI = kron(XYI,XB[nalleles-1])
		YYVI = kron(YYI,XB[nalleles-1])

		VI[XYVI] += rcorr
		VI[YYVI] -= rcorr

		VA[XYVA] += rcorr
		VA[YYVA] -= rcorr

		V[XYVA] += rcorr
		V[YYVA] -= rcorr

		for i = 1:nalleles-2

			XYVA = kron(BeyeA[i],kron(XYG,XB[nalleles-1-i])) | kron(XBA[i],kron(XYG,Beye[nalleles-1-i]))
			YYVA = kron(BeyeA[i],kron(YYG,XB[nalleles-1-i])) | kron(XBA[i],kron(YYG,Beye[nalleles-1-i]))
			XYVI = kron(BeyeI[i],kron(XYG,XB[nalleles-1-i])) | kron(XBI[i],kron(XYG,Beye[nalleles-1-i]))
			YYVI = kron(BeyeI[i],kron(YYG,XB[nalleles-1-i])) | kron(XBI[i],kron(YYG,Beye[nalleles-1-i]))

			VI[XYVI] += rcorr
			VI[YYVI] -= rcorr

			VA[XYVA] += rcorr
			VA[YYVA] -= rcorr

			V[XYVA] += rcorr
			V[YYVA] -= rcorr
		end

		XYVA = kron(XBA[nalleles-1],XYG)
		YYVA = kron(XBA[nalleles-1],YYG)
		XYVI = kron(XBI[nalleles-1],XYG)
		YYVI = kron(XBI[nalleles-1],YYG)

		VI[XYVI] += rcorr
		VI[YYVI] -= rcorr

		VA[XYVA] += rcorr
		VA[YYVA] -= rcorr

		V[XYVA] += rcorr
		V[YYVA] -= rcorr

		# Compute initial condition for Inactive state distribution
		# Steady state distribution is null space of V
		pss = nullspace(V)
		pss /= sum(pss)
		# println(pss)

		# Marginalize to obtain steady state distribution of 1 allele and nalleles-1 alleles
		nAnm1 = nG^(nalleles-1)
		pssnm1 = Array{Float64,1}(nAnm1)

		for i = 1:nAnm1
			pssnm1[i] = sum(pss[(i-1)*nA+1:i*nA])
		end
		# println(pssnm1)

		y = pss[kron(maskzb0z1,fill(true,nG^(nalleles-1)))]

		SIinit[kron(maskIinit,fill(true,nG^(nalleles-1)))] = y
		SIinit /= sum(SIinit)

		maskI = kron(maskI,ones(Int,nG^(nalleles-1)))

		SAinit = kron(SAinit,pssnm1)
		maskA = kron(maskA,ones(Int,nG^(nalleles-1)))

	else
		VA = TA
		VI = TI

		pss = nullspace(T)
		pss = pss/sum(pss)

		SIinit = [pss[mask]/sum(pss[mask]);zeros(div(length(SIinit),2))]

	end

	# Compute survival probability

	# Declare survival and dwell time distributions
	ntime = length(t)
	SIj = zeros(ntime)

	TIeiv = eigfact(VI)
	TIvects = TIeiv[:vectors]
	TIvals = TIeiv[:values]

	TIparams = TIvects\SIinit

	for j in find(maskI)
		for i = 1:size(TIvects,2)
			SIj += real(TIparams[i]*TIvects[j,i]*exp(TIvals[i]*t))
		end
	end

	SIj /= SIj[end]

	return SIj

end
