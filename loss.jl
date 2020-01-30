

function transition(n::Int,zeta::Int,r::Vector{Float64},nhist::Int)

	nT = (n+1)*2^zeta

	# Declare Transition matrices
	T = zeros(nT,nT)
	B = zeros(nT,nT)
	GT = zeros(nT,nT)
	GB = zeros(nT,nT)

	# Transition rates
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

	# Generate T = A + B transition matrix and B transition matrix
	for i=1:n+1, z=1:2^zeta, ip=1:n+1, w=1:2^zeta
		a = i + (n+1)*(z-1)
		b = ip + (n+1)*(w-1)
		zdigits = digits(z-1,base=2,pad=zeta)
		wdigits = digits(w-1,base=2,pad=zeta)
		z1 = zdigits[1]
		w1 = wdigits[1]
		zzeta = zdigits[zeta]
		wzeta = wdigits[zeta]
		zbar1 = zdigits[2:zeta]
		wbar1 = wdigits[2:zeta]
		zbarzeta = zdigits[1:zeta-1]
		wbarzeta = wdigits[1:zeta-1]
		T[a,b] =  (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)) + nu[zeta+1]*float((i==ip)*(zbarzeta==wbarzeta)*((zzeta==0)-(zzeta==1))*(wzeta==1))
		B[a,b] =  nu[zeta+1]*float((i==ip)*(zbarzeta==wbarzeta)*(zzeta==0)*(wzeta==1))
		GT[a,b] =  (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)) + nu[zeta+1]*float((i==ip)*(zbarzeta==wbarzeta)*((zzeta==0)-(zzeta==1))*(wzeta==1))
		GB[a,b] =  nu[zeta+1]*float((i==ip)*(zbarzeta==wbarzeta)*(zzeta==0)*(wzeta==1))
		for j = 1:zeta-1
			zbarj = zdigits[[1:j-1;j+2:zeta]]
			wbarj = wdigits[[1:j-1;j+2:zeta]]
			zj = zdigits[j]
			zj1 = zdigits[j+1]
			wj = wdigits[j]
			wj1 = wdigits[j+1]
			T[a,b] += nu[j+1]*float((i==ip)*(zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
		end
	end

	return T, B, nT, nu[end]

end

function mtransition(T::Matrix,B::Matrix,nT::Int,total::Int)

	S = zeros(total,total)
	Sminus = similar(S)
	Splus = similar(S)

	# Generate matrices for m transitions
	for m=1:total,mp=1:total
		S[m,mp] = -(mp-1)*float(m==mp)
		Sminus[m,mp] = float(m==mp+1)
		Splus[m,mp] = (mp-1)*float(m==mp-1)
	end

	M = kron(S,Matrix{Float64}(I,nT,nT)) + kron(Matrix{Float64}(I,total,total),T-B) + kron(Sminus,B) + kron(Splus,Matrix{Float64}(I,nT,nT))
	M[(total-1)*nT+1:end,(total-1)*nT+1:end] .+= B

	return M
end


function ss(n::Int,zeta::Int,rin::Vector{Float64},nhist::Int)
	# Steady state mRNA distribution of generalized telegraph model

	r = copy(rin)

	reactiondecay = 2*n + 2*zeta + 2
	r ./= rin[reactiondecay]

	T,B, nT, nu = transition(n,zeta,r,nhist)
	# A = T - B

	M = mtransition(T,B,nT,nhist)

	p = nullspace(M)[:,1]

	p /= sum(p)

	mhist = zeros(nhist)

	# Marginalize over G states
	for n in 1:nhist
		i = (n-1)*nT
		mhist[n] = sum(p[i+1:i+nT])
	end

	mhist = abs.(mhist)

	return mhist' * log.(mhist)

end

loss(x) = ss(1,1,x,10)
