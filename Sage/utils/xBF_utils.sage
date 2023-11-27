
def wt(x):
	return len([i for i in x if i != 0])

# ref: https://enac.hal.science/hal-02614017/document (pg.8)
def xBF(H,s,n,k,w,w_e,tau,N):
	t = s
	x = zero_vector(GF(2),n)
	e = zero_vector(GF(2),n-k)
	round = 0 
	while True:
		y = zero_vector(GF(2),n)
		for i in range(n):
			count = 0
			for j in range(n-k):
				if (t[j] == 1 and H[j,i]==1):
					count = count + 1
				if count >= tau:
					#print("count: ", count)
					y[i] = 1
		x = x + y
		#print(x)
		t = t + y * (H.transpose())
		if (wt(t) < w_e or round > N):
			break
		round += 1
	if (round <= N):
		return (True,x,s-x*(H.transpose()))
	# OSS: non è detto che il peso della soluzione sia w. (pseudocodice sbagliato)
	else: 
		return (False,x,s-x*(H.transpose()))
	


def test():
	k = 20
	n = 50
	MS = MatrixSpace(GF(2), k,n, sparse=True)
	M =  MS.random_element(density=0.1)
	G = block_matrix([[1,M]]) # generator matrix of the code
	H = block_matrix([[M.transpose(),1]]) # parity-check matrix of the code
	randomMessage = random_vector(GF(2),k)
	randomCodeword = randomMessage * G
	#creo un errore di qualche bit
	VS = MatrixSpace(GF(2), 1,n, sparse=True)
	v =  VS.random_element(density=0.1)
	receivedCodeord = randomCodeword + v

	res, = xBF(H,s,n,k,w,w_e,tau,N)