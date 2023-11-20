# description: Lee&Brickell complexity
# input: dimension of H (n,r), weight (w)
# output: complexity of L&B ISD
def LeeBrickelComplexity(n,r,w):
	k = n - r
	best_LB = 100000
	NC_w = (binomial(n,w)) / (2^r + 0.) # codeword estimate
	for p in range(30):
		p_ISD = binomial(k, p) * binomial(n-k, w-p) / binomial(n,w)
		p_star_ISD = NC_w * p_ISD
		if p_star_ISD>1:
			p_star_ISD=1
		t_solve = binomial(k,p)
		T = (n^3 + t_solve) / (p_star_ISD*1.) # complexity
		T = log(T*1.)/log(2.)
		if (T < best_LB):
			best_LB = T
			best_p = p
	return (best_LB, best_p)


# description: Stern Complexity
# input: dimension of H (n,r), weight (w)
# output: complexity of Stern ISD
def SternComplexity(n,r,w):
	k = n - r 
	best_S = 100000
	NC_w2 = binomial(n,w) / (2^r + 0.)
	for ell in range(30):
		for p in range(20):
			p_nkw = binomial(round(k+ell)/2,p)^2*binomial(n-k-ell,w-2*p)/binomial(n,w);
			p_star2 = NC_w2 * p_nkw
			if p_star2>1:
				p_star2=1
			L = binomial(round((k+ell)/2),p)
			T2 = (n^3 + L+L^2/(2^ell)) / (p_star2*1.)
			T2 = log(T2*1.)/log(2.)
			if T2 < best_S:
				best_S = T2
				best_p = p
				best_ell = ell
	return (best_S, best_p, best_ell)