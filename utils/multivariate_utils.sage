# description: characteristic function of x
# input: x
# output: {0,1}
def ch(x):
    if x != 0:
        return 1
    else:
        return 0



def createSystem(R,n,m):
	S = []
	for i in range(m):
		# Take random non-zero element r
		r = R._random_nonzero_element(randrange(n+1), randrange(n+ n*(n-1)/2 +1))
		#r = R._random_nonzero_element(2, randrange(n+ n*(n-1)/2 +1))
		# In the following lines: set all nonzero coeff to 1
		# (the problem is that r could be of the form x0^2*x1, which does not make sense in F2)
		if (r.degree() > 0):
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
		else: 
			i = i-1
	return S



def quadr(R,S,m):
	# We compute quadr(S)
	S1 = [[] for i in range(m)]	#Will contain the new (equiv) system, with all
								#polynomials lowered down to deg 2, in m different
								#lists (each containing the rearrangement of the
								#i-th polynomial)
	for i in range(m):
		poly = S[i]
		newPoly = 0
		for monomial in poly.monomials():
			if (monomial.degree() >= 2):
				nonzeroIndex = [index for index, value in enumerate(monomial.exponents()[0]) if value != 0]
				d = len(nonzeroIndex)
				for k in range(d-1):	#We need d-2 new variables
					R = PolynomialRing(GF(2), 'x', R.ngens()+1)
					if (k == 0):
						f_tmp = R.gens()[-1] + R.gens()[nonzeroIndex[0]]*R.gens()[nonzeroIndex[1]]
					else:
						f_tmp = R.gens()[-1] + R.gens()[-2]*R.gens()[nonzeroIndex[k+1]]
					S1[i].append(f_tmp)
				newPoly += R.gens()[-1] #* R.gens()[nonzeroIndex[-1]]
			else:
				newPoly += monomial 
		S1[i].append(newPoly)
	return R, S1



def removeOverlap(R, S1, m):
	S2 = [[] for i in range(m)]
	V = set()					#Set of variables appearing S1
	for i in range(m):	#Cycle over the m groups of polynomials
		S2[i] = S1[i][:-1]	#Copy all the polynomials in S1 of the form 'xy + z'
		# Cycle over the polynomials of the form 'xy + z'
		for j in range(len(S1[i][:-1])):
			var = S1[i][j].variables()
			for x in var:
				# Check whether a var was already seen, eventually substituting it
				if(x in V):
					R = PolynomialRing(GF(2),'x', R.ngens()+1)
					x_1 = R.gens()[-1]
					# Substitute x with x_1
					if(var == S2[i][j].variables()):
						S2[i][j] = S2[i][j].subs({x : x_1})	
					else:
						temp = S2[i][j].variables()[0]			#temp == x
						#print(f'Test: {temp==x}')
						S2[i][j] = S2[i][j].subs({temp : x_1})
					# Add the new relation
					S2[i].append(x + x_1)
				else:
					V.add(x)
		S2[i].append(S1[i][-1])	#Add the linear polynomial
	return R, S2



def linearFix(R, S2, m):
	S3 = [[] for i in range(m)]
	for i in range(m):	#Cycle over the m groups of polynomials
		S3[i] = S2[i][:-1]	#Copy all the polynomials of the form 'xy + z' and
							#'x + y'
		poly = S2[i][-1]	#Consider the other linear polynomial
		# Check whether poly is already in the desired form 'x + y + z + delta'
		var = poly.variables()	#Variables in poly
		d = len(var)			#NB: d = number of non-constant monomials
		if(d <= 3):
			S3[i].append(poly)
		else:
			for k in range(d-2):	#We need d-2 new variables
				if(k == 0):
					R = PolynomialRing(GF(2), 'x', R.ngens()+1)
					f = var[0] + var[1] + R.gens()[-1]			#= x_0 + x_1 + y_0
				elif(k == d-3):
					f = R.gens()[-1] + var[d-2] + var[d-1]		#= y_{d-3} + x_{d-2} + x_{d-1}
				else:
					R = PolynomialRing(GF(2), 'x', R.ngens()+1)
					f = R.gens()[-2] + var[k+1] + R.gens()[-1]	#= y_{k-1} + x_{k+1} + y_k
				S3[i].append(f)
	return R, S3


def printSystem(S):
	[print(p, end=' = 0\n') for p in S]
	print("*"*50)