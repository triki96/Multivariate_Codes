reset()
from random import randrange
from traceback import print_tb
from sage.all import *
load_attach_path('./utils/')
load('ISD_utils.sage')
load('multivariate_utils.sage')



#We create a random multivariate system
debug = False
if debug:
	m = 1 # number of equations
	n = 5 # number of variables
	set_random_seed(3)
else:
	m = 3 # number of equations
	n = 10 # number of variables

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

R, S3 = linearFix(R,S2,m)

print("New system (linear fix): \n")
for i in range(m):
	tmp = [print(p, end=' = 0\n') for p in S3[i]]
print("*"*50)


##################################################################################
#------------------------------Complexity Estimates------------------------------#
##################################################################################

#We count the number of equation we ended up with

q = 0 	#Number of quadratic equations
l = 0 	#Number of linear equations

for i in range(m):
	temp = len([p for p in S3[i] if p.degree() > 1])
	
	q += temp
	l += len(S3[i]) - temp

n = 10 * q
r = 7 * q - l
w = 3 * q

best_LB, best_p = LeeBrickelComplexity(n,r,w)
print("Complexity with Lee&Brickell ( p =",best_p,"): ", ceil(best_LB), "security bit")
print("*"*50)

best_S, best_p, best_ell = SternComplexity(n,r,w)
print("Complexity with Stern: ( p =",best_p,", l = ",best_ell,"):", ceil(best_S), "security bit")
print("*"*50)
