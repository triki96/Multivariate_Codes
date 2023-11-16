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
debug = False

if debug:
	m = 2 # number of equations
	n = 5 # number of variables

	set_random_seed(3)
else:
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
print("Original system: \n")
[print(p, end=' = 0\n') for p in S]
print("*"*50)


# We compute quadr(S)
S1 = [[] for i in range(m)]	#Will contain the new (equiv) system, with all
							#polynomials lowered down to deg 2, in m different
							#lists (each containing the rearrangement of the
							#i-th polynomial)

for i in range(m):
	poly = S[i]
	newPoly = 0

	for monomial in poly.monomials():
		if (monomial.degree() > 2):
			nonzeroIndex = [index for index, value in enumerate(monomial.exponents()[0]) if value != 0]
			d = len(nonzeroIndex)
		
			for k in range(d-2):	#We need d-2 new variables
				R = PolynomialRing(GF(2), 'x', R.ngens()+1)
		
				if (k == 0):
					f_tmp = R.gens()[-1] + R.gens()[nonzeroIndex[0]]*R.gens()[nonzeroIndex[1]]
		
				else:
					f_tmp = R.gens()[-1] + R.gens()[-2]*R.gens()[nonzeroIndex[k+1]]
		
				S1[i].append(f_tmp)
		
			newPoly += R.gens()[-1] * R.gens()[nonzeroIndex[-1]]
		
		else:
			newPoly += monomial 
	
	S1[i].append(newPoly)

print("New system (deg 2): \n")
for i in range(m):
	[print(p, end=' = 0\n') for p in S1[i]]
print("*"*50)


# Now we have to focus on polynomials of degree 2 which are not of the form
# xy + z = 0
S2 = [[] for i in range(m)]

for i in range(m):	#Cycle over the m groups of polynomials
	S2[i] = S1[i][:-1]	#Copy all the polynomials in S1 of the form 'xy + z'
	
	poly = S1[i][-1]	#Get the last one of each group, which by contruction
						#is the only one (possibly) not of the form 'xy + z' 

	if(poly.degree() == 1):
		S2[i].append(poly)

	else:
		newPoly = 0

		for mono in poly.monomials():
			if(mono.degree() == 2):
				R = PolynomialRing(GF(2),'x', R.ngens()+1)

				x = mono.variables()[0]
				y = mono.variables()[1]
				z = R.gens()[-1]

				S2[i].append(x*y + z)
				newPoly += z

			else:
				newPoly += mono

		S2[i].append(newPoly)

	#print(f'S2 = {S2}')

print("New system (xy+z = 0): \n")
for i in range(m):	
	[print(p, end=' = 0\n') for p in S2[i]]
print("*"*50)


# Removing overlapping variables in equations of the form "xy + z"
S3 = [[] for i in range(m)]
V = set()					#Set of variables appearing S2

for i in range(m):	#Cycle over the m groups of polynomials
	S3[i] = S2[i][:-1]	#Copy all the polynomials in S2 of the form 'xy + z'

	# Cycle over the polynomials of the form 'xy + z'
	for j in range(len(S2[i][:-1])):
		var = S2[i][j].variables()

		for x in var:
			# Check whether a var was already seen, eventually substituting it

			if(x in V):
				R = PolynomialRing(GF(2),'x', R.ngens()+1)

				x_1 = R.gens()[-1]
				
				# Substitute x with x_1
				if(var == S3[i][j].variables()):
					S3[i][j] = S3[i][j].subs({x : x_1})	
				else:
					temp = S3[i][j].variables()[0]			#temp == x
					#print(f'Test: {temp==x}')
					S3[i][j] = S3[i][j].subs({temp : x_1})

				# Add the new relation
				S3[i].append(x + x_1)

			else:
				V.add(x)

	S3[i].append(S2[i][-1])	#Add the linear polynomial

print("New system (non-overlapping): \n")
for i in range(m):
	[print(p, end=' = 0\n') for p in S3[i]]
print("*"*50)


# Now we work on linear polynomials such as f := x0 + x2 + x9 + 1 and we
# transform them according to the MPS paper

S4 = [[] for i in range(m)]

for i in range(m):	#Cycle over the m groups of polynomials
	S4[i] = S3[i][:-1]	#Copy all the polynomials of the form 'xy + z' and
						#'x + y'

	poly = S3[i][-1]	#Consider the other linear polynomial

	# Check whether poly is already in the desired form 'x + y + z + delta'
	var = poly.variables()	#Variables in poly
	d = len(var)			#NB: d = number of non-constant monomials
	
	if(d <= 3):
		S4[i].append(poly)
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
			
			S4[i].append(f)

print("New system (linear fix): \n")
for i in range(m):
	[print(p, end=' = 0\n') for p in S4[i]]
print("*"*50)


##################################################################################
#------------------------------Complexity Estimates------------------------------#
##################################################################################

#We count the number of equation we ended up with

q = 0 	#Number of quadratic equations
l = 0 	#Number of linear equations

for i in range(m):
	temp = len([p for p in S4[i] if p.degree() > 1])
	
	q += temp
	l += len(S4[i]) - temp

n = 10 * q
r = 7 * q - l
w = 3 * q

best_LB, best_p = LeeBrickelComplexity(n,r,w)
print("Complexity with Lee&Brickell ( p =",best_p,"): ", ceil(best_LB), "security bit")
print("*"*50)

best_S, best_p, best_ell = SternComplexity(n,r,w)
print("Complexity with Stern: ( p =",best_p,", l = ",best_ell,"):", ceil(best_S), "security bit")
print("*"*50)