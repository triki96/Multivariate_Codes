reset()
from random import randrange
from traceback import print_tb
from sage.all import *
import numpy as np

load_attach_path('./utils/')
load('ISD_utils.sage')
load('multivariate_utils.sage')
load('LB_utils.sage') # needed for Stern ISD
load('list_sorting_utils.sage') # needed for Stern ISD
load('stern_utils.sage') # needed for ISD
load('xBF_utils.sage')


#We create a random multivariate system
debug = False
if debug:
	m = 1 # number of equations
	n = 5 # number of variables
	set_random_seed(3)
else:
	m = 2 # number of equations
	n = 5 # number of variables

R = PolynomialRing(GF(2),'x', n) #Polynomial Ring

S = createSystem(R,n,m)

print("*"*50)
print("Original system: \n")
[print(p, end=' = 0\n') for p in S]
print("*"*50)

R, S1 = quadr(R,S,m)

print("New system (deg 2): \n")
for i in range(m):
	tmp = [print(p, end=' = 0\n') for p in S1[i]]
print("*"*50)

# Removing overlapping variables in equations of the form "xy + z"
R, S2 = removeOverlap(R, S1, m)

print("New system (non-overlapping): \n")
for i in range(m):
	tmp = [print(p, end=' = 0\n') for p in S2[i]]
print("*"*50)

# Now we work on linear polynomials such as f := x0 + x2 + x9 + 1 and we
# transform them according to the MPS paper

# R, S3 = linearFix(R,S2,m)
S3 = S2

print("New system (linear fix): \n")
for i in range(m):
	tmp = [print(p, end=' = 0\n') for p in S3[i]]
print("*"*50)


##################################################################################
#------------------------------Creating The Matrix-------------------------------#
##################################################################################

#We count the number of equation we ended up with

q = 0 	#Number of quadratic equations
l = 0 	#Number of linear equations

for i in range(m):
	temp = len([p for p in S3[i] if p.degree() > 1])
	q += temp
	l += len(S3[i]) - temp

#G_tmp = matrix([[1,0,0,1,1,0,0,1,1,1],[0,1,0,0,0,1,1,1,1,1],[0,0,1,1,1,1,1,1,1,1]])
H_tmp = matrix([[1,0,1],[1,0,1],[0,1,1],[0,1,1],[1,1,1],[1,1,1],[1,1,1]])
H_tmp = block_matrix([[H_tmp, 1]])

H = H_tmp
for i in range(q-1):
	H = block_diagonal_matrix(H, H_tmp)

#calcolo la prima parte della sindrome
syndrome = []
epsilon = vector(GF(2), [0,0,0,0,0,0,0,1,1,1])
for i in range(q):
	syndrome = vector(np.append(syndrome,epsilon))
syndrome = H * syndrome

M_tmp = matrix([[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0]])
M = M_tmp
for i in range(q-1):
	M = block_diagonal_matrix(M, M_tmp)

for i in range(m):
	temp = [p for p in S3[i] if p.degree() == 1] # prendo solo i polinomi lineari
	for f in temp:
		# creo un vettorel lungo n con 1 sui coeff giusti e lo aggiungo alla matrice
		polyToAdd = sum(vector(v) for v in f.exponents())
		padLength = len(R.gens()) - len(polyToAdd) # padding fino a H.nrows
		polyToAdd = vector(np.pad(polyToAdd, (0, padLength), 'constant'))
		polyToAdd = polyToAdd * M
		H = H.insert_row(H.nrows(),polyToAdd)
		syndrome = vector(np.append(syndrome, f.constant_coefficient()))

# H è la matrice di parità sulla quale vogliamo fare ISD
		

##################################################################################
#------------------------------Complexity Estimates------------------------------#
##################################################################################


n = 10 * q
r = 7 * q + l
w = 3 * q

#print('rank check: ',H.rank(), '=?=', r, "\n") # controllo se il rango è diminuito


best_LB, best_p = LeeBrickelComplexity(n,r,w)
print("Complexity with Lee&Brickell ( p =",best_p,"): ", ceil(best_LB), "security bit")
print("*"*50)

best_S, best_p, best_ell = SternComplexity(n,r,w)
print("Complexity with Stern: ( p =",best_p,", l = ",best_ell,"):", ceil(best_S), "security bit")
print("*"*50)


##################################################################################
#-----------------------Compute & Transform Solutions----------------------------#
##################################################################################

check, ISD_sol, ISD_syndrome = xBF(H,syndrome,n,n-r,w,0,3,200)

if check:
	sol = [i[1] for i in enumerate(ISD_sol) if ((i[0]%10==0 or i[0]%10==1 or i[0]%10==2) and (i[0] < 7 * q))]
	print("Solution to multivariate system: ", sol)
else:
	print("No solution found")