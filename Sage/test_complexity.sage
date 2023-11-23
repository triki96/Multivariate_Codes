# Creo un sistema random, faccio numTests test e vedo la complessità cosa viene fuori
# Creo un sistema UOV-style, faccio numTests test e vedo la complessità cosa viene fuori
reset()
from random import randrange
from traceback import print_tb
from sage.all import *
import numpy as np
from cysignals.alarm import alarm, AlarmInterrupt, cancel_alarm

load_attach_path('./utils/')
load('ISD_utils.sage')
load('multivariate_utils.sage')
load('LB_utils.sage') # needed for Stern ISD
load('list_sorting_utils.sage') # needed for Stern ISD
load('stern_utils.sage') # needed for ISD


#We create a random multivariate system
m = 3 # number of equations
n = 10 # number of variables
v = 3
o = n - v
complexityCounterRandom = zero_vector(256)
complexityCounterUOV = zero_vector(256)
numTests = 1001

for i in range(numTests):
	#print("new try: ", i)
	try:
		alarm(60) # due minuti di tempo massimo
		R = PolynomialRing(GF(2),'x', n) #Polynomial Ring
		S = createSystem(R,n,m) # Create random system
		R, S1 = quadr(R,S,m) # Put system in quadr form
		R, S2 = removeOverlap(R, S1, m) # remove overlapping variables
		S3 = S2 # We don't fix the linear part
		H,q = createMatrixCode(S3) # Create The Matrix
		#-Complexity Estimates-#
		nn = H.ncols()
		rr = H.rank() #r = 7 * q + l
		ww = 3 * q

		best_S, best_p, best_ell = SternComplexity(nn,rr,ww)
		best_S = ceil(best_S)
		if (best_S > 255):
			best_S = 255
		complexityCounterRandom[best_S] += 1
		#print(i,"-randomCounter: ", complexityCounterRandom)
		if (i % 100 == 0):
			print(i,"-randomCounter: ", complexityCounterRandom)
	except (AlarmInterrupt):
		#print("Didnt manage")
		cancel_alarm()
		i = i-1 # retry with new system
	else:
		#print("OK")
		cancel_alarm()






for i in range(numTests):
	try:
		alarm(60) # due minuti di tempo massimo
		R = PolynomialRing(GF(2),'x', n) #Polynomial Ring
		S = generateUOVPK(2,m,n,v,o)
		R, S1 = quadr(R,S,m) # Put system in quadr form
		R, S2 = removeOverlap(R, S1, m) # remove overlapping variables
		S3 = S2 # We don't fix the linear part
		H,q = createMatrixCode(S3) # Create The Matrix

		#Complexity Estimates#
		nn = H.ncols()
		rr = H.rank() #r = 7 * q + l
		ww = 3 * q

		best_S, best_p, best_ell = SternComplexity(nn,rr,ww)
		best_S = ceil(best_S)
		if (best_S > 255):
			best_S = 255
		complexityCounterUOV[best_S] += 1
		if (i % 100 == 0):
			print(i,"-randomCounter: ", complexityCounterUOV)
	except (AlarmInterrupt):
		cancel_alarm()
		i = i-1 # retry with new system
	else:
		cancel_alarm()
