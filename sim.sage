reset()
from random import randrange
from traceback import print_tb
from sage.all import *


# description: characteristic function of x
# input: x
# output: {0,1}
def ch(x):
    if x != 0:
        return 1
    else:
        return 0
	

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
	

#We create a random multivariate system
max_degree = 2
m = 3 # number of equations
n = 10 # number of variables
R = PolynomialRing(GF(2),'x', n) #Polynomial Ring
S = []
for i in range(m):
	# Take random non-zero element r
	r = R._random_nonzero_element(randrange(n+1), randrange(n+ n*(n-1)/2 +1))
	#r = R._random_nonzero_element(2, randrange(n+ n*(n-1)/2 +1))
	# In the following lines: set all nonzero coeff to 1
	# (the problem is that r could be of the form x0^2*x1, which does not make sense in F2)
	r_toAppend = 0
	for j in range(len(r.monomials())):
		monomial = r.monomials()[j]
		exponents = zero_vector(n)
		for k in range(n):
			if (monomial.exponents()[0][k] > 0):
				exponents[k] = 1
		newMonomial = prod(var**(ch(exp)) for var, exp in zip(R.gens(), exponents))
		r_toAppend += newMonomial
	S.append(r_toAppend)


print("*"*50)
print("Original system: \n", S)
print("*"*50)


# We compute quadr(S)
S1 = [] #will contain the new (equiv) system, with all polynomials lowered down to deg 2
for i in range(m):
	poly = S[i]
	newPoly = 0
	for j in range(len(poly.monomials())):
		monomial = poly.monomials()[j]
		if (monomial.degree() > 2):
			nonzeroIndex = [index for index, value in enumerate(monomial.exponents()[0]) if value != 0]
			d = len(nonzeroIndex)
			for k in range(d-2):
				# we need d-2 new variables
				R = PolynomialRing(GF(2),'x', R.ngens()+1)
				if (k == 0):
					f_tmp = R.gens()[-1] + R.gens()[nonzeroIndex[k]]*R.gens()[nonzeroIndex[k+1]]
				else:
					f_tmp = R.gens()[-1] + R.gens()[-2]*R.gens()[nonzeroIndex[k+1]]
				S1.append(f_tmp)
			newPoly += R.gens()[-1] * R.gens()[nonzeroIndex[-1]]
		else:
			newPoly += monomial 
	S1.append(newPoly)


print("New system (deg 2): \n", S1)
print("*"*50)


#Now we have to focus on polynomials of degree 2 which are not of the form
# xy + z = 0
S2 = []
for i in range(len(S1)): # we cycle over all polynomials
	newPoly = 0
	if (S1[i].degree() != 1):
		monomials = S1[i].monomials()
		numMonomials = len(monomials)
		if (numMonomials > 2):
			for j in range(numMonomials): #cycle over all monomials
				if (monomials[j].degree() == 2):  # we just focus on deg 2
					R = PolynomialRing(GF(2),'x', R.ngens()+1)
					f_tmp = R.gens()[-1] + (monomials[j].variables())[0] * (monomials[j].variables())[1]
					S2.append(f_tmp)
					newPoly += R.gens()[-1]
				else:
					newPoly += monomials[j]
		else:
			newPoly += S1[i]
	else:
		newPoly += S1[i]
	S2.append(newPoly)


print("New system: (xy+z = 0): \n", S2)
print("*"*50)


#Now we work on linear polynomials such as f := x0 + x2 + x9 + 1 and we transform
#them according to the MPS paper
S3 = []
for i in range(len(S2)): # cycle over all polynomials
	f_tmp = S2[i] # this will be added to S3
	if (S2[i].degree() != 1): # we focus on deg 1 polynomials
		S3.append(f_tmp)
	else:
		numMonomials =  len(S2[i].exponents())
		if (numMonomials <= 3): # in particular we focus on linear polynomials with more that 3 monomials
			S3.append(f_tmp)
		else:
			if (S2[i].constant_coefficient() == 1): #see if the polynomial contains +1
				#print(S[i], " - caso in cui il polinomio contiene il termine noto")
				if (numMonomials > 4):
					variables = S2[i].variables()
					d = len(variables)
					for k in range(d-2):
					# we need d-2 new variables
						R = PolynomialRing(GF(2),'x', R.ngens()+1)
						if (k == 0):
							f_tmp = R.gens()[-1] + variables[0] + variables[1]
						elif (k == d-3):
							f_tmp = R.gens()[-1] + variables[k+1] + 1
						else:
							f_tmp = R.gens()[-1] + R.gens()[-2] + variables[k+1]
						S3.append(f_tmp)
				else:
					S3.append(f_tmp)
			else:
				#print(S2[i], " - caso in cui il polinomio NON contiene il termine noto")
				variables = S2[i].variables()
				d = len(variables)
				for k in range(d-2):
				# we need d-2 new variables
					if (k == 0):
						R = PolynomialRing(GF(2),'x', R.ngens()+1)
						f_tmp = R.gens()[-1] + variables[0] + variables[1]
					elif (k == d-3):
						f_tmp = R.gens()[-1] + variables[k+1] + variables[k+2]
					else:
						R = PolynomialRing(GF(2),'x', R.ngens()+1)
						f_tmp = R.gens()[-1] + R.gens()[-2] + variables[k+1]
					S3.append(f_tmp)
				

print("New system (linear fix): \n", S3)
print("*"*50)


##################################################################################
#------------------------------Complexity Estimates------------------------------#
##################################################################################


#We count the number of equation we ended up with
q = len([i for i in S3 if i.degree() > 1]) #number of quadratic equations
l = len(S3) - q #number of linear equations
n = 10 * q
r = 7 * q - l
w = 3 * q

best_LB, best_p = LeeBrickelComplexity(n,r,w)
print("Complexity with Lee&Brickell ( p =",best_p,"): ", ceil(best_LB), "security bit")
print("*"*50)

best_S, best_p, best_ell = SternComplexity(n,r,w)
print("Complexity with Stern: ( p =",best_p,", l = ",best_ell,"):", ceil(best_S), "security bit")
print("*"*50)