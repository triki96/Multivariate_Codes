reset()
from random import randrange
from traceback import print_tb
from sage.all import *

def generateUOVPK(q,m,n,v,o):
	R = PolynomialRing(GF(q),'x', n) #Polynomial Ring
	
	# S is the private UOV key
	F = []
	for k in range(m):
		# Take random non-zero element r
		r = 0
		for i in range(v):
			for j in range(v):
				coeff = randrange(q)
				r = r + coeff * R.gens()[i] * R.gens()[j]
		for i in range(v):
			for j in range(o):
				coeff = randrange(q)
				r = r + coeff * R.gens()[i] * R.gens()[v+j]
		F.append(r)

	# T is the set of affine trasformations
	T = []
	for k in range(n):
		r = R._random_nonzero_element(1, randrange(n+ n*(n-1)/2 +1))
		T.append(r)

	#Masked system
	P = []
	for k in range(m):
		P.append(F[k](tuple(T[i] for  i in range(n))))

	return P


#######################################################################################


def test():
	q = 2
	m = 3 		# number of equations
	n = 10	 	# number of variables
	v = 5			# number of vinegard variables
	o = n - v	# number of oil variables

	P =generateUOVPK(q,m,n,v,o)
	print(P)